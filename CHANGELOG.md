# Changelog

## (unreleased)

### Other

* Update README.md.


## 0.0.26 (2024-07-04)

### Other

* Merge pull request #25 from TomHarrop/0.0.26.

  add quality trimming to final step

* Add quality trimming to final step.


## 0.0.25 (2024-04-18)

### Other

* Check sample names for double underscores.

* Update graphs.


## 0.0.24 (2023-10-16)

### Other

* Have to tell clumpify to use tmp.


## 0.0.23 (2023-10-16)

### Other

* Switch to clumpify for duplicate counting.


## 0.0.22 (2023-10-15)

### Other

* Release 0.0.22.

* Don't try to scale resources (issue #22)


## 0.0.21 (2023-10-15)

### Other

* Scale cores for duplicates.

* Use quay API for latest tag.


## 0.0.20 (2023-10-13)

### Other

* Actually learn the learnings.


## 0.0.19 (2023-10-13)

### Other

* TIL python exponents.


## 0.0.18 (2023-10-12)

### Other

* Merge pull request #21 from TomHarrop/dedupe_again.

  don't count containments, it takes too long

* Don't count containments, it takes too long.


## 0.0.17 (2023-10-09)

### Other

* Merge pull request #20 from TomHarrop/fix_mem.

  Set mem and threads to global resources to avoid the oom killer

* Rst path.

* Also set threads to be global.

* Linting.

* Readme.

* Fix mem_mb param.

* Update graphs.


## 0.0.16 (2023-10-06)

### Other

* Merge pull request #19 from TomHarrop/stats.

  Stats

* Linting.

* Also plot duplication stats.

* Add stats scripts.

* Keep full cutadapt stats.

* Stats_typo.

* Added cutadapt stat parsing.

* Start parsing json.

* Lint.

* Cutadapt stats in json format.


## 0.0.15 (2023-09-14)

### Other

* Merge pull request #18 from TomHarrop/0.0.15.

  hook into snakemake logger

* Hook into snakemake logger.


## 0.0.14 (2023-09-13)

### Other

* Merge pull request #17 from TomHarrop/0.0.14.

  use periods for delimiting read number so that we can parse the sampl…

* Use periods for delimiting read number so that we can parse the sample name.

* Don't need to print tags.

* Env.new_release.

* Update builder.yml.

* Check release logic.

* Merge pull request #16 from TomHarrop/TomHarrop-patch-1.

  try using python only

* Update readme.

* Update README.md.

* Just create release, no more docker stuff.

* Try using python only.


## 0.0.13 (2023-09-08)

### Other

* Merge pull request #15 from TomHarrop/0.0.13.

  add test for missing files (fixes #14)

* Remove comment.

* Add test for missing files (fixes #14)

* Tidy up filenams.


## 0.0.12 (2023-09-01)

### Other

* Merge pull request #13 from TomHarrop/0.0.12.

  Added test for special characters in filenames (fixes #11)

* Added test for special characters in filenames (fixes #11)


## 0.0.11 (2023-09-01)

### Other

* Merge pull request #12 from TomHarrop/0.0.11.

  memory for repair

* Memory for repair.


## 0.0.10 (2023-09-01)

### Other

* Merge pull request #10 from TomHarrop/0.0.10.

  Allow parallel demux but be stricter with RAM (#7)

* Allow parallel demux but be stricter with RAM (#7)


## 0.0.9 (2023-09-01)

### Other

* Merge pull request #9 from TomHarrop/0.0.9.

  Hard code RAM and prevent more than one demux job (#7)

* Hard code RAM and prevent more than one demux job (#7)


## 0.0.8 (2023-08-31)

### Other

* Merge pull request #8 from TomHarrop/0.0.8.

  Try setting threads and using less ram for barcode check (issue #7)

* Try setting threads and using less ram for barcode check (issue #7)

* Readme typos.


## 0.0.7 (2023-08-17)

### Other

* Merge pull request #5 from TomHarrop/0.0.7.

  0.0.7

* Add thread advice.

* Simpler threading, hopefully higher throughput.


## 0.0.6 (2023-08-16)

### Other

* Merge pull request #4 from TomHarrop/0.0.6.

  0.0.6

* Try to improve parallelism.

* Try to improve parallelism.

* Update readme.

* Update readme.

* Update readme.

* Update readme.

* Update readme.

* Update readme.

* Update readme.

* Update readme.


## 0.0.5 (2023-08-11)

### Other

* Dont lock workdir.


## 0.0.4 (2023-08-11)

### Other

* Include script.


## 0.0.3 (2023-08-11)

### Other

* Merge pull request #3 from TomHarrop/0.0.3.

  update thread logic

* Update thread logic.


## 0.0.2 (2023-08-11)

### Other

* Merge pull request #2 from TomHarrop/0.0.2.

  try installing to site-packages

* Try installing to site-packages.

* Try installing to site-packages.

* Merge pull request #1 from TomHarrop/package.

  Initial release

* Tidy.

* Tidy.

* Tidy.


## 0.0.1 (2023-08-11)

### Other

* Tag.

* Tidy.

* Manual version.

* Rebuild.

* VERSION.

* Initial package build attempt.

* Update docker actions.

* Update release-tagging.

* Get planned release tag for test.

* Don't build the Singularity container.

* Token needs permission?

* Remove deprecated action.

* Update builder.yml.

* Initial commit.

* Update builder.yml.

* Force in changes from template.

* Need sudo to rm android and ghc.

* Handle free space better.

* Still need CR_PAT for pull???

* Token.

* Layer caching.

* Only build IF not latest tag.

* Update template Dockerfile.

* Try release.

* Variable imagename.

* Specify file for version to work.

* Test with version.

* Needs to be one workflow to keep variables.

* Try getting lowercase tags.

* Try getting lowercase tags.

* Try getting tags.

* Try getting tags.

* Try getting tags.

* Try getting tags.

* Actually build masurca.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Test login.

* Try pull after push:

* Rerun docker test.

* Test time.

* Build branch.

* Initial commit.

* Create LICENSE.

* Add setup.py.

* Tempdir.

* Readme.

* Modify stats location.

* Test all samples.

* Add external pipeline.

* Add mask and trim.

* Move functions.

* Update design in readme.

* Add logic for external barcode only samples.

* Write barcorde file to disk.

* Working first demux.

* Working first demux.

* Semi-working demux.

* Semi-working demux.

* Initial demux step to check barcodes.

* Parse metadata.

* Update strategy.

* Start design.

* Start design.

* Initial commit.


