# ExoFOP-TESS submission package: TIC 71316387 (CTOI 71316387.01)

*Prepared 2026-10-16. **Nothing has been posted.** This is a draft for the account holder to review and post. AI-run;
no professional astronomer has reviewed it.*

The package has two parts:

1. **An observing note** in ExoFOP house style.
2. **An 8-page supporting report (PDF)** to upload under the same data tag.

Both state that the signal is an eclipsing binary at twice the CTOI period, most likely an M-dwarf pair bound to the
G star. They rest on the refereed campaign in [`../README.md`](../README.md).

## Files

| File | What it is |
|---|---|
| [`observing_note.txt`](observing_note.txt) | The note, plain ASCII, 1,222 characters, one line break |
| [`observing_note_bulk.txt`](observing_note_bulk.txt) | The same note with `<br>` for the line break, for ExoFOP's bulk-upload route |
| [`TIC71316387.01O-xx20261016-EB-analysis.pdf`](TIC71316387.01O-xx20261016-EB-analysis.pdf) | Supporting report: 8 pages, 3 figures, 5 tables |
| [`build/`](build/) | Scripts that make the figures, the eclipse-time table and the PDF, with their outputs |
| [`review/review_log.md`](review/review_log.md) | Every reviewer finding and what was done about it |
| [`review/review_findings.json`](review/review_findings.json) | Raw reviewer and drafting output |

## The note

```
CTOI 71316387.01 is an EB at 2x the CTOI period, not a planet: P = 1.30367228 +/- 0.00000019 d, primary T0 = BJD_TDB 2460165.67042 +/- 0.00018 (TESS-SPOC S1/28/68, TESScut S95/105/106); the CTOI epoch is a secondary eclipse. Depths 5770/2600 ppm, deeper on the same parity in all six sectors. The G star is not in the EB: its Gaia DR3 RV is constant (15 RVS transits over 946 d, chi2 p = 0.28), vs K of 18 km/s or more for any stellar EB partner. Colour-corrected difference-image centroids put the source within 1.24 arcsec of the G star (95%); all Gaia neighbours (18 arcsec and beyond) excluded. The TESS centroid colour shift and SuperWASP/NEOWISE W1 depths show a red (early-M, possibly late-K) eclipsed star; the eclipse shape with MS relations gives an ~M2V + ~M4V pair with ~2.6% of TESS flux (not dynamical). Likely a bound triple (chance alignment disfavoured ~5000-14000x for baseline priors, ~200x worst case), unconfirmed. Suggest FP (blended EB).
To confirm: AO/speckle (predicted delta-Ks ~2.6, delta-m(832 nm) ~3.8 mag; within 0.9 arcsec, possibly within 0.1 arcsec), then CPM or M-dwarf RVs. No G-star PRV needed. AI-assisted analysis (Claude Code); full report: TIC71316387.01O-xx20261016-EB-analysis.pdf
```

**Why it looks like this.** The note follows a study of all 33,156 public ExoFOP observing notes (2026-10-16):

- **Length.**
  - The median note is 215 characters; the 90th percentile is 632.
  - The closest analogues, multi-part EB and source-attribution notes by independent vetters, run 647–2,226 characters
    (median 916, IQR 804–1,375).
  - Notes over ~1,800 characters often needed correction notes later.
- **Form.** The verdict comes first and the period is stated against the listed one. The new T0 is labelled BJD_TDB,
  each line of evidence gets one clause, and the disposition is suggested rather than decreed.
- **ASCII only.** Arcsec is written out because ExoFOP's export drops `"`. "Below/within" replaces `<`, and `+/-`
  replaces `±`, because non-Latin-1 characters become `¿`.
- **Firsts.** This would be the first observing note on this target and on any RAVEN CTOI. No existing note mentions AI
  assistance. The one-clause disclosure is a choice, not a convention; remove it if you prefer, since the PDF discloses in
  full.

## Before posting: placeholders to fill

| Where | Placeholder | Replace with |
|---|---|---|
| PDF file name and the last words of the note | `xx` in `TIC71316387.01O-xx20261016-EB-analysis.pdf` | Your initials. ExoFOP's pattern is `<target>O-<initials><UT date>-<info>.pdf`. Use the upload date if it is not 2026-10-16 |
| PDF title page and metadata | `[Author name] (ExoFOP-TESS user [username])` | Your name and ExoFOP username |

To rebuild the PDF with your details (system python3 with reportlab; the figures are already in `build/figs/`):

```
cd tess-search/exofop_leads/TIC71316387_campaign/exofop_submission/build
REPORT_AUTHOR="Your Name" REPORT_AUTHOR_LINE="Your Name (ExoFOP-TESS user yourname)" \
REPORT_NAME=TIC71316387.01O-AB20261016-EB-analysis.pdf python3 build_report.py   # writes the PDF into build/
```

Or ask Claude to do it.

## How to post (from ExoFOP's own documentation)

1. **Log in** to ExoFOP-TESS.
2. **Pick one named data tag** for the note and the file, for example `20261016_<username>_EBvetting_71316387`.
   - The description part may use letters, numbers and hyphens only.
   - Named tags shared by a note and its files are the norm; auto-tags almost never carry files.
3. **Upload the PDF** on the TIC 71316387 page: Files, then "Add new". Settings:
   - Type: O (Other).
   - Tag: the one above.
   - Group: blank. Proprietary period: 0 (only TFOP members can set one).
   - Suggested description: "Supporting report: CTOI 71316387.01 is an EB at 2x the CTOI period, most likely a bound
     M-dwarf pair; AI-assisted analysis of public data".
4. **Add the observing note** on the same target with the same tag.
   - Type it into the editor or paste it as plain text, not from a chat window: pasted formatting is stored verbatim.
   - Enter the single line break before "To confirm:" with the Enter key.
5. **Optional courtesy.** The ExoFOP conduct policy asks users to contact uploaders before using their unpublished data.
   We used only public data and the published RAVEN paper. Still, telling the CTOI submitter (`lafargamagro`, RAVEN team)
   would be courteous.

Do not create a new CTOI. ExoFOP's candidate rules say a period that is an integer multiple of an existing candidate's
belongs in a note.

## How it was checked

All of the following steps were done by AI agents.

**The note.**
- Three independent drafts in different house styles were each scored by three judges.
- A synthesis followed, then a fact-check of every number against the campaign results and the raw Gaia and CTOI data.
- A second, sceptic's fact-check could not run because the account's usage limit was reached. It was done by hand
  instead.

**The report.** Three reviewers went over it, yielding about 75 findings:
- an audit that recomputed every number;
- an adversarial TFOP-style referee;
- a presentation review.

Every finding was either applied or answered in [`review/review_log.md`](review/review_log.md). The planned reference
check and per-finding verification could not run for the same reason. They were done by hand: each change was re-checked
against the source files, and the RAVEN citation was confirmed by web search.

**Changes the review produced.**
- R★ is now 1.08 R☉ (Gaia FLAME), not 1.07.
- Fig. 3a now plots the field-star calibration that the quoted colour actually comes from.
- The V-band depth prediction is now a range.
- The W1 range for a G/F star is wider.
- Data products are named correctly (TESS-SPOC FFI).
- A Gaia parallax zero-point note was added.
- A secondary-eclipse argument now explains why the G star cannot be eclipsed.

**Repo fix.** The review also found that two campaign result files had been swapped. They are fixed in the repo.

## Limits

- The bound-triple conclusion is probabilistic: no resolved image, common proper motion or M-dwarf RVs exist yet.
- The masses, light fraction and spectral types rest on main-sequence relations.
- Section 8 of the report lists the open issues. They include a 1.6–1.9σ temperature tension, shallow-eclipse colours
  that disagree with the model, and an unexplained faint VHS source at 1.1″.
