********************
[0.2.0] - XXXX-XX-XX
********************

**New features**:

**Breaking changes**:

- Bumped SampleData file format version to 4.0 then to 5.0, as a result of adding flags
  on individuals and setting tskit.UNKNOWN_TIME for sites where freq is used as time.
- The ancestors tree sequence now contains the real alleles and not
  0/1 values as before.
- Times for undated sites now use frequencies (0..1), not as counts (1..num_samples),
  and are now stored as tskit.UNKNOWN_TIME, then calculated on the fly when ancestors
  are generated.
- The SampleData file no longer accepts the ``inference`` argument to add_site.
  This functionality has been replaced by the ``exclude_positions`` argument
  to the ``infer`` and ``generate_ancestors`` functions.

**Bugfixes**:

- Individuals and populations in the sampledata file are kept in the returned tree
  sequence, even if they are not referenced by any sample. The individual and population
  ids are therefore guaranteed to stay the same between the sample data file and the
  inferred tree sequence. (:pr:`348`)

********************
[0.1.5] - 2019-09-25
********************

**Breaking changes**:

- Bumped SampleData file format version to 2.0, then to 3.0 as a result of the additions
  below. Older SampleData files will not be readable and must be regenerated.

- Users can specify variant ages, via ``sample_data.add_sites(... , time=user_time)``.
  If not ``None``, this overrides the default time position of an ancestor, otherwise
  ancestors are ordered in time by using the frequency of the derived variant (#143).
  This addition bumped the file format to 2.0

- Change "age" to "time" to match tskit/msprime notation, and to avoid confusion
  with the age since birth of an individual (#149). Together with the 2 changes below,
  this addition bumped the file format to 3.0.

- Add the ability to record user-specified times for individuals, and therefore
  the samples contained in them (currently ignored during inference). Times are
  added using ``sample_data.add_individual(... , time=user_time)`` (#190). Together
  with the changes above and below, this addition bumped the file format to 3.0.

- Change ``tsinfer.UNKNOWN_ALLELE`` to ``tskit.MISSING_DATA`` for marking unknown regions
  of ancestral haplotypes (#188) . This also involves changing the allele storage to a
  signed int from ``np.uint8`` which matches the tskit v0.2 format for allele storage
  (see https://github.com/tskit-dev/tskit/issues/144). Together with the 2 changes above,
  this addition bumped the file format to 3.0.

**New features**:

- Map non-inference sites onto the tree by using the built-in tskit
  ``map_mutatations`` method. With further work, this should allow triallelic sites
  to be mapped (#185)

- The default tree sequence returned after inference when ``simplify=True`` retains
  unary nodes (i.e. simplify is done with ``keep_unary=True``. This tends to result
  in better compression.


********************
[0.1.4] - 2018-12-12
********************

Bugfix release.

- Fix issue caused by upstream changes in numcodecs (#136).

********************
[0.1.3] - 2018-11-02
********************

Release corresponding to code used in the preprint.

********************
[0.1.2] - 2018-06-18
********************

Minor update to take advantage of msprime 0.6.0's Population and Individual
objects and fix various bugs.


**Breaking changes**:

- Bumped SampleData file format version to 1.0 because of the addition
  of individuals and populations. Older SampleData files will not be
  readable and must be regenerated.

- Changed the order of the ``alleles`` and ``genotypes`` arguments to
  SampleData.add_site.

**New features**:

- Sample and individual metadata now handled correctly.

- Added --no-simplify option to CLI and simplify=True option to infer function.

- Better handling of missing files (raises correct exceptions).

- tsinfer list now presents basic information for .trees files.

**Bug fixes**:

- Degenerate examples with zero inference sites are now rooted (#44)

- Simplify=False results in tree sequence with correct sample nodes.
