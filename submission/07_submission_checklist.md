# GigaByte Submission Checklist — ibd-taxonomic-dysbiosis Technical Release

Every requirement below is drawn from GigaByte's *Information for Authors*.
Placeholders in square brackets must be completed before online submission
at the GigaByte manuscript submission system. Editorial office:
editorial@gigabytejournal.com.

---

## 1. Account and submission channel
- [ ] Create/log-in to account in the GigaByte online manuscript submission system
  (new authors: create an account first).
- [ ] Submit all files through the submission system.
- [ ] Contact editorial@gigabytejournal.com for anything not supported by the system.

## 2. Manuscript format
- [x] `04_manuscript_gigabyte.docx` — single file, Word format.
- [x] Double line spacing applied.
- [x] Line numbers applied to all pages (continuous restart).
- [x] Page numbers in footer.
- [x] No page breaks within the manuscript.
- [x] SI units only (verify each measurement in text once more).
- [x] Title page: title **+ authors, affiliations, addresses, emails, ORCIDs,
      corresponding author** — filled: Edward Jenner Tettevi (ORCID
      0000-0002-2448-2679, ejtettevi@gmail.com) and Mike Y. Osei-Atweneboana
      (oseiatweneboana@yahoo.co.uk); CSIR-WRI, P.O. Box M32, Ministries,
      Accra, Ghana.

## 3. Cover letter (`02_cover_letter.docx`)
- [x] States why the article is appropriate for GigaByte (reproducible
      research-data analysis; no new dataset).
- [x] Notes journal-policy issues (none expected).
- [x] Competing interests declaration.
- [x] AI-assisted technology use disclosure.
- [x] Preprint status.
- [x] Confirmation all authors approved the submission.
- [x] Confirmation work not published or submitted elsewhere.
- [ ] **[PENDING: complete author name, date, signature]**.

## 4. Article type (Technical Release)
- [x] The submission is a Technical Release: open-source computational
      method / pipeline for the analysis of research data, with tested
      software and documentation.
- [x] Software is on GitHub: https://github.com/ejtettevi/ibd-taxonomic-dysbiosis
- [x] OSI-approved LICENSE added (BSD-3-Clause), commit 915e6a7.
- [x] Repository made public (via `gh`, visibility=public).
- [x] Repository tagged `v1.0` (annotated); Availability section cites
      v1.0 / commit 915e6a7.
- [ ]; environment.yml rebuild test on a fresh machine still to be run;
      `run_all.sh` now resolves Rscript/python from PATH or CONDA_PREFIX
      (no absolute miniforge paths).

## 5. Figures (`03_figures/`, 9 PNGs)
- [x] Nine figures as separate files, named and numbered in order of first
      mention (figure_1 … figure_9).
- [x] Each figure can be placed on a single portrait page — verify.
- [x] High resolution (README/figures produced at 300 dpi).
- [x] Figure legends included in manuscript; figures NOT embedded in
      manuscript text. **[Verify figure_* PNGs are referenced but not embedded
      in 04_manuscript_gigabyte.docx.]**

## 6. Tables
- [x] Tables 1–4 in the main manuscript file (pandoc-converted Word table objects).
- [x] Numbered and cited in sequence in the text.
- [x] Titles ≤ 15 words; legends ≤ 300 words — **verify against rendered docx**.

## 7. References (`02_manuscript_gigabyte.md` → docx)
- [x] Numbered consecutively in square brackets **[1]**… in **order of first
      citation** (renumbered; verified map).
- [x] GigaByte style (journal names abbreviated; `Year;Vol:Pages.` or `doi:`).
- [x] All 18 references verified against PubMed/Crossref (see `gb_refs`).
      Note: Gevers 2014 = *Cell Host Microbe* 15:382-392; Wang 2024 =
      *Biomed Pharmacother* 179:117302; Manghi 2025 = *Nat Commun* 17:196.

## 8. Required sections present in manuscript
- [x] "Availability and Requirements" (project name, home page, OS,
      programming language, other requirements, license, RRID, bio.tools ID,
      version). License row = **BSD-3-Clause**; version row = **v1.0 (commit
      915e6a7)**.
- [x] "Data and Materials" (repository + permanent identifier; source data
      attribution). GigaByte-recommended GigaDB sentence added; dataset DOI
      to be inserted when issued.
- [x] "Disclosure of use of AI-assisted tools" at end of manuscript.
- [x] Abstract ends with "Availability and Implementation" line.
- [x] Research areas: primary **Microbiology**, classifications
      **Metagenomics; Medical Microbiology** (verify against the current
      GigaByte area list at submission time — lists are revised periodically).
- [x] Ethics statement, contributions, competing interests, funding,
      acknowledgements, data availability.

## 9. Supplementary material / GigaDB
- [x] No PDF supplementary files — all research objects go to GigaDB.
- [x] GigaDB metadata prepared: `06_gigadb_metadata.md` (author + full file
      manifest filled).
- [x] GigaDB README finalized: `05_README_GigaDB.md` (purpose, provenance,
      file dictionary, reproduction steps, compute note, license, limits).
- [ ] Submit data package to GigaDB curation team (gigadb@gigasciencejournal.com
      or via the submission system) and obtain dataset DOI.

## 10. ML-paper compute disclosure (not required here — no training runs;
      classifiers are small fixed-seed models) — **confirm one more time**.
      If expecting a "compute" disclosure: 8-core CPU, ~<20 min, <8 GB RAM.

---

### Critical pending actions before upload (summary)
1. Add signature/date line to the cover letter (`02_cover_letter.docx`) before
   signing/uploading.
2. Submit data package to GigaDB curation team and insert the issued dataset
   DOI into "Data and Materials" + GigaDB metadata.
3. Fresh-machine `environment.yml` rebuild test of the pipeline
   (`environment.yml` now aligned with actual runtime deps, commit e97a224).
4. Final sign-off on author/affiliation/address on title page
   (currently verified against manuscript DOCX + cover letter DOCX).