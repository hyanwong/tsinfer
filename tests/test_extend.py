import numpy as np
import pytest

import tsinfer


def assert_variants_equal(vars1, vars2):
    assert vars1.num_sites == vars2.num_sites
    assert vars1.num_samples == vars2.num_samples
    for var1, var2 in zip(vars1.variants(), vars2.variants()):
        assert var1.alleles == var2.alleles
        assert np.all(var1.genotypes == var2.genotypes)


class TestExtend:
    @pytest.mark.parametrize("num_samples", range(1, 5))
    @pytest.mark.parametrize("num_sites", range(1, 5))
    def test_single_binary_haplotype_one_generation(self, num_samples, num_sites):
        with tsinfer.SampleData(sequence_length=num_sites) as sd:
            for j in range(num_sites):
                sd.add_site(j, [1] * num_samples)
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend(np.arange(num_samples))
        assert_variants_equal(ts, sd)

    @pytest.mark.parametrize("k", range(1, 5))
    def test_single_binary_haplotype_k_generations(self, k):
        num_sites = 5
        num_samples = 4
        with tsinfer.SampleData(sequence_length=num_sites) as sd:
            for j in range(num_sites):
                sd.add_site(j, [1] * (num_samples * k))

        extender = tsinfer.SequentialExtender(sd)
        for _ in range(k):
            ts = extender.extend(np.arange(num_samples) * k)
        assert_variants_equal(ts, sd)

    def test_single_haplotype_4_alleles(self):
        num_sites = 3
        with tsinfer.SampleData(sequence_length=num_sites) as sd:
            for j in range(num_sites):
                sd.add_site(j, [0, 1, 2, 3], alleles="ACGT")

        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend(np.arange(4))
        assert_variants_equal(ts, sd)

    @pytest.mark.parametrize("k", range(4, 9))
    def test_single_site_4_alleles_rotating(self, k):
        genotypes = np.zeros(k, dtype=int)
        for j in range(k):
            genotypes[j] = j % 4
        with tsinfer.SampleData(sequence_length=1) as sd:
            sd.add_site(0, genotypes, alleles="ACGT")

        extender = tsinfer.SequentialExtender(sd)
        for j in range(k):
            ts = extender.extend([j])
        assert ts.num_mutations == 3
        assert_variants_equal(ts, sd)

    @pytest.mark.parametrize("num_generations", range(1, 5))
    @pytest.mark.parametrize("samples_per_generation", [1, 2, 13])
    @pytest.mark.parametrize("num_sites", [1, 4, 10, 100])
    def test_random_data(self, num_generations, samples_per_generation, num_sites):
        rng = np.random.default_rng(42)
        num_samples = num_generations * samples_per_generation
        with tsinfer.SampleData(sequence_length=num_sites) as sd:
            for j in range(num_samples):
                sd.add_individual(ploidy=1, metadata={"ind_id": j})
            for j in range(num_sites):
                genotypes = rng.integers(0, 4, size=num_samples)
                sd.add_site(j, genotypes, alleles="ACGT")

        extender = tsinfer.SequentialExtender(sd)
        offset = 0
        for _ in range(num_generations):
            next_offset = offset + samples_per_generation
            ts = extender.extend(np.arange(offset, next_offset))
            assert ts.num_samples == next_offset
            offset = next_offset

        assert ts.num_sites == sd.num_sites
        assert ts.num_samples == sd.num_samples
        for var1, var2 in zip(ts.variants(alleles=("A", "C", "G", "T")), sd.variants()):
            assert var1.alleles == var2.alleles
            assert np.all(var1.genotypes == var2.genotypes)
        for j, u in enumerate(ts.samples()):
            assert ts.node(u).metadata == {"ind_id": j}

    def test_single_sample_metadata(self):
        with tsinfer.SampleData(sequence_length=1) as sd:
            sd.add_individual(ploidy=1, metadata={"x": 1})
            sd.add_site(0, [1])
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0])
        assert_variants_equal(ts, sd)
        assert ts.node(ts.samples()[0]).metadata == {"x": 1}

    @pytest.mark.parametrize("num_generations", range(1, 5))
    def test_stick(self, num_generations):
        # We have a stick tree where the single mutation for a given site
        # happens on one branch and they accumulate over time.
        H = np.zeros((num_generations, num_generations), dtype=int)
        a = np.zeros(num_generations, dtype=int)
        for j in range(num_generations):
            a[j] = 1
            H[j] = a
        with tsinfer.SampleData(sequence_length=num_generations) as sd:
            for j in range(num_generations):
                sd.add_site(j, H[:, j])
        extender = tsinfer.SequentialExtender(sd)
        for j in range(num_generations):
            ts = extender.extend([j])
            assert ts.num_samples == j + 1
        assert ts.num_mutations == num_generations
        assert ts.num_edges == num_generations + 1

    @pytest.mark.parametrize("num_generations", range(1, 5))
    def test_all_zeros(self, num_generations):
        # all the haplotypes are 0s and should just copy directly from
        # the same root.
        a = np.zeros(2 * num_generations, dtype=int)
        with tsinfer.SampleData(sequence_length=num_generations) as sd:
            sd.add_site(0, a)
        extender = tsinfer.SequentialExtender(sd)
        for j in range(num_generations):
            ts = extender.extend([2 * j, 2 * j + 1])
            # assert ts.num_samples == 2 * j + 1
        assert ts.num_mutations == 0
        assert ts.num_trees == 1
        tree = ts.first()
        parents = {tree.parent(u) for u in ts.samples()}
        assert len(parents) == 1

    @pytest.mark.parametrize("num_generations", range(1, 5))
    def test_provenance(self, num_generations):
        a = np.zeros(num_generations, dtype=int)
        with tsinfer.SampleData(sequence_length=num_generations) as sd:
            sd.add_site(0, a)
        extender = tsinfer.SequentialExtender(sd)
        for j in range(num_generations):
            ts = extender.extend([j])
            assert ts.num_provenances == 1


class TestExtendPathCompression:
    def example(self):
        with tsinfer.SampleData(sequence_length=4) as sd:
            sd.add_site(0, [0, 1, 1, 1])
            sd.add_site(1, [0, 1, 1, 1])
            sd.add_site(2, [1, 0, 1, 1])
            sd.add_site(3, [1, 0, 2, 1], alleles=("0", "1", "2"))
            return sd

    def test_simple_path_compression_case(self):
        sd = self.example()
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0, 1])
        # NOTE we'd really like to get rid of this vestigial node 0 but
        # the low-level code won't work without it, so until it's
        # gone it's simplest to just live with it and update the test
        # cases later.

        # 2.00┊  0  ┊
        #     ┊  ┃  ┊
        # 1.00┊  1  ┊
        #     ┊ ┏┻┓ ┊
        # 0.00┊ 2 3 ┊
        #     0     4
        assert ts.num_trees == 1
        assert ts.num_nodes == 4
        assert ts.first().parent_dict == {2: 1, 3: 1, 1: 0}

        ts = extender.extend([2, 3])
        # 3.00┊   0   ┊   0   ┊
        #     ┊   ┃   ┊   ┃   ┊
        # 2.00┊   1   ┊   1   ┊
        #     ┊ ┏━┻┓  ┊ ┏━┻┓  ┊
        # 1.00┊ 2  3  ┊ 3  2  ┊
        #     ┊    ┃  ┊    ┃  ┊
        # 1.00┊    6  ┊    6  ┊
        #     ┊   ┏┻┓ ┊   ┏┻┓ ┊
        # 0.00┊   4 5 ┊   4 5 ┊
        #     0       2       4

        assert ts.num_trees == 2
        assert ts.num_nodes == 7
        assert ts.node(6).flags == tsinfer.NODE_IS_PC_ANCESTOR
        assert ts.first().parent_dict == {2: 1, 3: 1, 1: 0, 6: 3, 4: 6, 5: 6}
        assert ts.last().parent_dict == {2: 1, 3: 1, 1: 0, 6: 2, 4: 6, 5: 6}
        assert_variants_equal(ts, sd)

    def test_simple_path_compression_case_no_pc(self):
        sd = self.example()

        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0, 1])
        assert ts.num_trees == 1
        assert ts.num_nodes == 4
        assert ts.first().parent_dict == {2: 1, 3: 1, 1: 0}

        ts = extender.extend([2, 3], path_compression=False)
        # 3.00┊   0   ┊   0   ┊
        #     ┊   ┃   ┊   ┃   ┊
        # 2.00┊   1   ┊   1   ┊
        #     ┊ ┏━┻┓  ┊ ┏━┻┓  ┊
        # 1.00┊ 2  3  ┊ 3  2  ┊
        #     ┊   ┏┻┓ ┊   ┏┻┓ ┊
        # 0.00┊   4 5 ┊   4 5 ┊
        #     0       2       4
        assert ts.num_trees == 2
        assert ts.num_nodes == 6
        assert ts.first().parent_dict == {2: 1, 3: 1, 1: 0, 4: 3, 5: 3}
        assert ts.last().parent_dict == {2: 1, 3: 1, 1: 0, 4: 2, 5: 2}
        assert_variants_equal(ts, sd)


class TestExtendIdenticalSequences:
    def test_single_site_one_generation(self):
        with tsinfer.SampleData(sequence_length=1) as sd:
            sd.add_site(0, [1, 1])
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0, 1])

        # 2.00┊  0  ┊
        #     ┊  ┃  ┊
        # 1.00┊  1  ┊
        #     ┊  ┃  ┊
        # 0.01┊  4  ┊
        #     ┊ ┏┻┓ ┊
        # 0.00┊ 2 3 ┊
        #     0     1
        assert ts.num_trees == 1
        assert ts.num_nodes == 5
        assert ts.first().parent_dict == {2: 4, 3: 4, 4: 1, 1: 0}
        assert ts.node(4).flags == tsinfer.NODE_IS_IDENTICAL_SAMPLE_ANCESTOR

        assert_variants_equal(ts, sd)

    def test_two_haplotypes_one_generation(self):
        alleles = ("A", "C", "G")
        with tsinfer.SampleData(sequence_length=2) as sd:
            sd.add_site(0, [1, 1, 2, 2], alleles=alleles)
            sd.add_site(1, [1, 1, 2, 2], alleles=alleles)
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0, 1, 2, 3])

        # 2.00┊    0    ┊
        #     ┊    ┃    ┊
        # 1.00┊    1    ┊
        #     ┊  ┏━┻━┓  ┊
        # 0.00┊  6   7  ┊
        #     ┊ ┏┻┓ ┏┻┓ ┊
        # 0.00┊ 2 3 4 5 ┊
        #     0         2
        assert ts.num_trees == 1
        assert ts.num_nodes == 8

        assert ts.first().parent_dict == {2: 6, 3: 6, 4: 7, 5: 7, 6: 1, 7: 1, 1: 0}
        assert ts.node(6).flags == tsinfer.NODE_IS_IDENTICAL_SAMPLE_ANCESTOR
        assert ts.node(7).flags == tsinfer.NODE_IS_IDENTICAL_SAMPLE_ANCESTOR

        assert_variants_equal(ts, sd)

    def test_two_haplotypes_two_generations(self):
        alleles = ("A", "C", "G")
        with tsinfer.SampleData(sequence_length=2) as sd:
            sd.add_site(0, [1, 1, 2, 2, 2, 2], alleles=alleles)
            sd.add_site(1, [1, 1, 2, 2, 2, 2], alleles=alleles)
        extender = tsinfer.SequentialExtender(sd)
        ts = extender.extend([0, 1, 2, 3])
        ts = extender.extend([4, 5])
        # We correctly see that there was a pre-existing exact match for
        # this haplotype and match against it.
        assert_variants_equal(ts, sd)
        # 3.00┊     0       ┊
        #     ┊     ┃       ┊
        # 2.00┊     1       ┊
        #     ┊  ┏━━┻━━┓    ┊
        # 1.00┊  6     7    ┊
        #     ┊ ┏┻┓ ┏━┳┻┳━┓ ┊
        # 1.00┊ 2 3 4 5 ┃ ┃ ┊
        #     ┊         ┃ ┃ ┊
        # 0.00┊         8 9 ┊
        #     0             2
        assert ts.first().parent_dict == {
            2: 6,
            3: 6,
            4: 7,
            5: 7,
            6: 1,
            7: 1,
            1: 0,
            8: 7,
            9: 7,
        }
        assert ts.node(6).flags == tsinfer.NODE_IS_IDENTICAL_SAMPLE_ANCESTOR
        assert ts.node(7).flags == tsinfer.NODE_IS_IDENTICAL_SAMPLE_ANCESTOR
