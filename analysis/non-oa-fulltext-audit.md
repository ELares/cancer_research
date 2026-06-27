# Non-OA full-text redistribution audit (#526)

Each `is_oa: false` corpus record with a PMCID, checked against the NCBI PMC Open Access subset service (`oa.fcgi`). Records NOT in the OA subset have their `## Full Text` stripped in place (abstract + metadata kept; the frozen corpus count is unchanged) via `scripts/strip_non_oa_fulltext.py --apply`.

- Non-OA records with a PMCID: **61**
- In the PMC OA subset (redistributable, KEPT): **27**
- NOT in the OA subset (copyrighted, STRIPPED): **34**

## Stripped (not in the OA subset)

| PMID | PMCID | oa_status | Journal |
|---|---|---|---|
| 19660870 | PMC2784186 | closed | Cancer treatment reviews |
| 20866097 | PMC2997921 | closed | Molecular pharmaceutics |
| 22033517 | PMC3280949 | closed | Science (New York, N.Y.) |
| 22542702 | PMC3413735 | closed | Bio Systems |
| 22586319 | PMC4211116 | closed | Cancer discovery |
| 22778154 | PMC3551628 | closed | Molecular cancer therapeutics |
| 23196890 | PMC3528107 | closed | Physical biology |
| 23749887 | PMC3709593 | unknown | Anticancer research |
| 24562770 | PMC4318538 | closed | Breast cancer research and treatment |
| 25003941 | PMC4137229 | closed | Physics in medicine and biology |
| 26079252 | PMC4808585 | closed | Critical reviews in clinical laboratory sciences |
| 26719576 | PMC4747838 | closed | Molecular cancer therapeutics |
| 26773162 | PMC4717912 | closed | Clinical cancer research : an official journal of the American Association for Cancer Research |
| 27693939 | PMC5108677 | closed | Bioelectrochemistry (Amsterdam, Netherlands) |
| 27815355 | PMC5413401 | closed | Clinical cancer research : an official journal of the American Association for Cancer Research |
| 28138868 | PMC5576028 | closed | Medical oncology (Northwood, London, England) |
| 28148839 | PMC5435119 | closed | Science translational medicine |
| 28493544 | PMC5647195 | closed | Journal of surgical oncology |
| 29144754 | PMC5821496 | closed | Nano letters |
| 29174271 | PMC5835831 | closed | Journal of shoulder and elbow surgery |
| 29180466 | PMC5811386 | closed | Cancer research |
| 29431697 | PMC5882515 | closed | Cancer discovery |
| 29492540 | PMC6340641 | closed | Acta neurochirurgica. Supplement |
| 29777637 | PMC6105452 | closed | Journal of oral pathology & medicine : official publication of the International Association of Oral Pathologists and the American Academy of Oral Pathology |
| 30240926 | PMC6289793 | closed | Photodiagnosis and photodynamic therapy |
| 30277116 | PMC6445778 | closed | Leukemia & lymphoma |
| 30339360 | PMC6702128 | closed | ACS applied materials & interfaces |
| 30341213 | PMC6279584 | closed | Cancer immunology research |
| 30415456 | PMC6414051 | closed | Journal of neuro-oncology |
| 30866031 | PMC7278092 | unknown | Oncology (Williston Park, N.Y.) |
| 30911535 | PMC6430577 | closed | Clinical and translational imaging |
| 31409607 | PMC6774877 | closed | Cancer immunology research |
| 31451760 | PMC6911768 | closed | Nature chemical biology |
| 31636445 | PMC7030949 | closed | Nature reviews. Nephrology |
