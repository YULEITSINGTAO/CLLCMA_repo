#!/usr/bin/env python3
"""Train DeepSurv (Cox-PH MLP) on fixed CV folds and write OOF risk scores.

Input TSV columns: PATIENT_ID, FOLD, TIME, EVENT, then numeric features.
Output TSV: PATIENT_ID, FOLD, RISK_SCORE
  RISK_SCORE is the network log-hazard (higher = higher risk), same direction
  as a Cox linear predictor.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import sys

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


META_COLS = {"PATIENT_ID", "FOLD", "TIME", "EVENT"}


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


class DeepSurvNet(nn.Module):
    """Katzman et al. DeepSurv: SELU MLP → scalar log-hazard."""

    def __init__(self, d_in: int, hidden: int = 32, dropout: float = 0.3):
        super().__init__()
        self.fc1 = nn.Linear(d_in, hidden)
        self.fc2 = nn.Linear(hidden, hidden)
        self.out = nn.Linear(hidden, 1, bias=False)
        self.dropout = nn.Dropout(dropout)
        for layer in (self.fc1, self.fc2, self.out):
            nn.init.normal_(layer.weight, mean=0.0, std=1.0 / np.sqrt(layer.weight.shape[1]))
            if layer.bias is not None:
                nn.init.zeros_(layer.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.dropout(F.selu(self.fc1(x)))
        x = self.dropout(F.selu(self.fc2(x)))
        return self.out(x).squeeze(-1)


def cox_ph_loss(log_h: torch.Tensor, time: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    """Negative mean Breslow partial log-likelihood. log_h higher = higher risk."""
    order = torch.argsort(time, descending=True)
    log_h = log_h[order]
    event = event[order]
    log_cumsum = torch.logcumsumexp(log_h, dim=0)
    n_events = event.sum()
    if n_events <= 0:
        return log_h.new_zeros(())
    return -((log_h - log_cumsum) * event).sum() / n_events


def read_table(path: str):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty table: {path}")
    feat_cols = [c for c in rows[0] if c not in META_COLS]
    pid = [r["PATIENT_ID"] for r in rows]
    fold = np.array([int(float(r["FOLD"])) for r in rows], dtype=np.int32)
    time = np.array([float(r["TIME"]) for r in rows], dtype=np.float32)
    event = np.array([int(float(r["EVENT"])) for r in rows], dtype=np.int32)
    x = np.array([[float(r[c]) for c in feat_cols] for r in rows], dtype=np.float32)
    x[~np.isfinite(x)] = 0.0
    return pid, fold, time, event, x


def zscore_fit(x: np.ndarray):
    mu = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    sd[~np.isfinite(sd) | (sd < 1e-8)] = 1.0
    return mu, sd


def zscore_apply(x: np.ndarray, mu: np.ndarray, sd: np.ndarray) -> np.ndarray:
    z = (x - mu) / sd
    z[~np.isfinite(z)] = 0.0
    return z


def event_stratified_val_idx(event: np.ndarray, frac: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    val = np.zeros(len(event), dtype=bool)
    for ev in (0, 1):
        idx = np.flatnonzero(event == ev)
        if idx.size == 0:
            continue
        rng.shuffle(idx)
        n_val = max(1, int(round(idx.size * frac)))
        n_val = min(n_val, idx.size - 1) if idx.size > 1 else 0
        if n_val > 0:
            val[idx[:n_val]] = True
    return val


def train_one_fold(
    x_train: np.ndarray,
    time_train: np.ndarray,
    event_train: np.ndarray,
    x_test: np.ndarray,
    seed: int,
    hidden: int,
    dropout: float,
    lr: float,
    weight_decay: float,
    max_epochs: int,
    patience: int,
    device: torch.device,
) -> np.ndarray:
    mu, sd = zscore_fit(x_train)
    x_tr = zscore_apply(x_train, mu, sd)
    x_te = zscore_apply(x_test, mu, sd)

    use_val = int(event_train.sum()) >= 10 and x_tr.shape[0] >= 40
    if use_val:
        val_mask = event_stratified_val_idx(event_train, 0.2, seed)
        if int(event_train[~val_mask].sum()) < 5 or int(event_train[val_mask].sum()) < 3:
            use_val = False

    if use_val:
        x_fit, t_fit, e_fit = x_tr[~val_mask], time_train[~val_mask], event_train[~val_mask]
        x_val, t_val, e_val = x_tr[val_mask], time_train[val_mask], event_train[val_mask]
    else:
        x_fit, t_fit, e_fit = x_tr, time_train, event_train
        x_val = t_val = e_val = None

    model = DeepSurvNet(x_fit.shape[1], hidden=hidden, dropout=dropout).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    def to_t(a, dtype):
        return torch.tensor(a, dtype=dtype, device=device)

    xt = to_t(x_fit, torch.float32)
    tt = to_t(t_fit, torch.float32)
    et = to_t(e_fit, torch.float32)
    if use_val:
        xv = to_t(x_val, torch.float32)
        tv = to_t(t_val, torch.float32)
        ev = to_t(e_val, torch.float32)

    best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    best_loss = float("inf")
    bad = 0
    epochs = max_epochs if use_val else min(max_epochs, 200)

    for _ in range(epochs):
        model.train()
        opt.zero_grad(set_to_none=True)
        loss = cox_ph_loss(model(xt), tt, et)
        if not torch.isfinite(loss):
            break
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()

        model.eval()
        with torch.no_grad():
            monitor = cox_ph_loss(model(xv), tv, ev) if use_val else loss
        monitor_v = float(monitor.detach().cpu())
        if np.isfinite(monitor_v) and monitor_v < best_loss - 1e-6:
            best_loss = monitor_v
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            bad = 0
        else:
            bad += 1
            if use_val and bad >= patience:
                break

    model.load_state_dict(best_state)
    model.to(device)
    model.eval()
    with torch.no_grad():
        risk = model(to_t(x_te, torch.float32)).detach().cpu().numpy().astype(np.float64)
    risk[~np.isfinite(risk)] = np.nan
    return risk


def run(input_path: str, output_path: str, seed: int, device_name: str) -> None:
    pid, fold, time, event, x = read_table(input_path)
    n = len(pid)
    risk = np.full(n, np.nan, dtype=np.float64)

    if device_name == "cpu":
        device = torch.device("cpu")
    elif device_name == "cuda" and torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    hidden = 32 if x.shape[1] >= 16 else 16
    folds = sorted(set(int(f) for f in fold))
    for f in folds:
        test = fold == f
        train = ~test
        if test.sum() < 5 or train.sum() < 30 or int(event[train].sum()) < 5:
            continue
        set_seed(seed + int(f))
        try:
            risk[test] = train_one_fold(
                x[train],
                time[train],
                event[train],
                x[test],
                seed=seed + int(f),
                hidden=hidden,
                dropout=0.3,
                lr=1e-3,
                weight_decay=1e-4,
                max_epochs=400,
                patience=40,
                device=device,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"Fold {f} failed: {exc}", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["PATIENT_ID", "FOLD", "RISK_SCORE"])
        for i in range(n):
            val = "" if not np.isfinite(risk[i]) else f"{risk[i]:.10g}"
            w.writerow([pid[i], int(fold[i]), val])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--device", default=os.environ.get("DEEPSURV_DEVICE", "cpu"))
    args = p.parse_args()
    run(args.input, args.output, args.seed, args.device)


if __name__ == "__main__":
    main()
