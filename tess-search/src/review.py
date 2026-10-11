"""Adversarial review of a candidate dossier by an independent model via OpenRouter.

Usage: review.py MODEL DOSSIER.md OUT.md [image.png ...]
Spending is checked against the OpenRouter key usage recorded before this task started; the call is
refused if this task would exceed its budget.
"""
import sys, json, base64, requests

BASELINE = 14.0151406   # key spend (USD) before this task, read 2026-10-11
CAP = 10.0              # budget for this task (USD)
URL = "https://openrouter.ai/api/v1"

SYSTEM = """You are an experienced exoplanet / eclipsing-binary vetting scientist reviewing work done by an AI agent.
Your job is to find reasons the claimed signal is NOT a new astrophysical discovery: instrumental systematics,
stellar variability, background or bound eclipsing binaries, contamination from neighbouring stars, period
aliases, look-elsewhere effects, catalogue matches that were missed, statistical errors, or flawed tests.
Be concrete and quantitative. Check the arithmetic. Say which claims are not supported by the evidence shown,
what additional test would most change your assessment, and give an overall verdict:
(a) likely real and new, (b) plausible but unproven, (c) likely false positive or already known."""


def spent():
    k = requests.get(f"{URL}/key", timeout=30).json()["data"]
    return float(k.get("usage", k.get("budgets", [{}])[0].get("spend_usd", 0.0)))


def review(model, text, images=(), max_tokens=24000):
    used = spent() - BASELINE
    if used > CAP - 0.5:
        raise SystemExit(f"budget nearly exhausted: ${used:.2f} of ${CAP:.2f} used")
    content = [{"type": "text", "text": text}]
    for p in images:
        b64 = base64.b64encode(open(p, "rb").read()).decode()
        content.append({"type": "image_url", "image_url": {"url": f"data:image/png;base64,{b64}"}})
    body = {"model": model, "max_tokens": max_tokens, "usage": {"include": True},
            "messages": [{"role": "system", "content": SYSTEM}, {"role": "user", "content": content}]}
    r = requests.post(f"{URL}/chat/completions", json=body, timeout=900)
    r.raise_for_status(); j = r.json()
    out = j["choices"][0]["message"]["content"]
    cost = j.get("usage", {}).get("cost")
    return out, cost, spent() - BASELINE


if __name__ == "__main__":
    model, dossier, out = sys.argv[1:4]
    text, cost, total = review(model, open(dossier).read(), sys.argv[4:])
    open(out, "w").write(f"# Review by {model}\n\n{text}\n")
    print(f"cost ${cost} ; task total so far ${total:.3f}")
