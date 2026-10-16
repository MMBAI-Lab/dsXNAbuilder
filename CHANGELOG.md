# Changelog

All notable changes to dsXNAbuilder. Versions follow [semantic versioning](https://semver.org/);
between releases the package version carries a `-dev` suffix. Every PDB file records the version
that built it, and the web app shows it in its footer.

## [Unreleased]

First public version.

### Builder
- DNA, RNA and DNA:RNA hybrid duplexes in 40 helical forms derived from the 3DNA fiber models
  (A, B, C, D, Z, B′, S and B*-DNA, A- and A′-RNA, hybrids), each nucleotide placed by the screw
  symmetry of its model.
- Both strands typed explicitly (5′→3′ over 3′→5′), with mismatches, overhangs and abasic sites.
- modXNA chemistry on either strand: 244 prebuilt residues (2′-deoxyribose, ribose, 2′-OMe, 2′-F,
  MOE and LNA sugars; A, C, G, T, U and 5-methylcytosine; phosphodiester and the two
  phosphorothioates; internal and terminal forms), transplanted onto the fiber scaffold with
  stereochemistry, bond and contact checks.
- Hydrogens from the Amber or modXNA residue templates, under the names tleap expects; flexible
  2′-O substituents turned out of contacts.
- Every finding reported with a stable code (errors stop the build, warnings do not).

### Output
- A self-describing PDB file (input, model, modified residues and data versions in its header),
  with residue names chosen by the user if wanted.
- Optional Amber/modXNA topology package (tleap input, libraries, frcmod, manifest with the
  command that rebuilds it).
- High-throughput lists in a FASTA-like format: up to 50 duplexes in the web app and 500 on the
  command line, returned as one zip with a report.

### Web app and command line
- Static web app: the Python package runs in the browser (Pyodide in a Web Worker); NGL viewer,
  modXNA library tab, reproducible links, documentation.
- Command-line tool `dsxnabuilder models | build` with the same results.

### Validation
- pytest suite, including base-pair and step parameters checked against CPPTRAJ `nastruct`;
  every library residue built with AmberTools (tleap, sander); end-to-end test of the web app.
