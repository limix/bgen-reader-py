def _make_genotypes(ploidy, nalleles):
    if ploidy > 0:
        tups = []
        for i in range(nalleles + 1):
            t = _make_genotypes(ploidy - 1, nalleles - i)
            for ti in t:
                nleft = nalleles - len(ti)
                tups += [nleft * [ploidy + 1] + ti]
        return tups
    else:
        return [nalleles * [ploidy + 1]]


def get_genotypes(ploidy, nalleles):
    t = _make_genotypes(ploidy - 1, nalleles)
    return [list(reversed(i)) for i in t]


def genotypes_to_allele_counts(genotypes):
    ploidy = genotypes[-1][0]
    counts = []
    for g in genotypes:
        count = [0] * ploidy
        for gi in g:
            count[gi - 1] += 1
        counts.append(count)
    return counts


# if __name__ == "__main__":
#     genotypes = _get_genotypes(3, 3)
#     for t in genotypes:
#         print(t)

#     print("-----------------")
#     counts = _convert_to_allele_counts(genotypes)
#     for t in counts:
#         print(t)

#     print("-------NEW----------")
#     genotypes = _get_genotypes(2, 3)
#     for t in genotypes:
#         print(t)

#     print("-----------------")
#     counts = _convert_to_allele_counts(genotypes)
#     for t in counts:
#         print(t)

#     print("-------NEW----------")
#     genotypes = _get_genotypes(2, 2)
#     for t in genotypes:
#         print(t)

#     print("-----------------")
#     counts = _convert_to_allele_counts(genotypes)
#     for t in counts:
#         print(t)
