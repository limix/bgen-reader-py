import datetime
import sys
import time
from contextlib import contextmanager


def get_genotypes(ploidy, nalleles):
    g = [_make_genotype(p, 1, nalleles) for p in ploidy]
    g = sorted([list(reversed(i)) for i in g])
    g = [list(reversed(i)) for i in g]
    return g


def genotypes_to_allele_counts(genotypes):
    nalleles = genotypes[-1][0]
    counts = []
    for g in genotypes:
        count = [0] * nalleles
        for gi in g:
            count[gi - 1] += 1
        counts.append(count)
    return counts


def _make_genotype(ploidy, start, end):
    tups = []
    if ploidy == 0:
        return tups
    if ploidy == 1:
        return [[i] for i in range(start, end + 1)]
    for i in range(start, end + 1):
        t = _make_genotype(ploidy - 1, i, end)
        for ti in t:
            tups += [[i] + ti]
    return tups


# Later: consider switching to tqdm package?
# best location?
@contextmanager
def _log_in_place(name, verbose, time_lambda=time.time, show_log_diffs=False):
    """
    Create a one-argument function to write messages to. All messages will be on the same line and
    will include time.

    """
    # LATER what if logging messages aren't suppose to go to stdout?
    t_wait = time_lambda()
    last_len = [0]  # We have to make this an array so that the value is by reference.
    last_message_hash = [None]
    every_printed = [False]  # Don't print the final newline if nothing is ever printed

    if not verbose:

        def do_nothing(message):
            pass

        yield do_nothing
    else:

        def writer(message):
            time_str = str(datetime.timedelta(seconds=time_lambda() - t_wait))
            if "." in time_str:
                time_str = time_str[
                    : time_str.index(".") + 3
                ]  # Time to the 1/100th of a sec
            s = "{0} -- time={1}, {2}".format(name, time_str, message)
            if show_log_diffs:
                message_hash = hash(message)
                if (
                    message_hash != last_message_hash[0]
                    and last_message_hash[0] is not None
                ):
                    sys.stdout.write("\n")
                last_message_hash[0] = message_hash
            # Pad with spaces to cover up previous message
            sys.stdout.write("{0}{1}\r".format(s, " " * max(0, last_len[0] - len(s))))
            sys.stdout.flush()
            every_printed[0] = True
            last_len[0] = len(s)

        yield writer

    if not verbose:
        return
    if every_printed[0]:
        sys.stdout.write("\n")
