# ftree

[![crates.io](https://img.shields.io/crates/v/ftree.svg)](https://crates.io/crates/ftree)
[![docs](https://docs.rs/ftree/badge.svg)](https://docs.rs/ftree)

A pure-rust(with zero dependencies, no-std) fenwick tree, for the efficient computation of dynamic prefix sums.

# Background

Did you ever have to keep track of a sum, and update it at the same time?

Let's say that you have an array that represents the lengths of some other containers:
[1, 6, 3, 9, 2]

What if you want to get the sum up until the n-th element? In the worst-case, this will take O(n) time. Updating on the other hand, is simply
a matter of incrementing at the specified index, at O(1).

A fenwick tree allows you to both get the sum and do updates in O(log n) time.

Moreover, let's assume that you want to get the index of the first value such that <= sum.

Without using a Fenwick tree, this would take (n * log n) time (a binary search with the sum being computed during each iteration). Using
one however, only takes O(log n) time. This might seem like a very niche need, and it is. It is utilized in the [indexset](https://crates.io/crates/indexset)
crate, a two-level B-Tree, to very efficiently support vector-like indexing by position.

# Performance

It's very performant. I have searched all over codeforces for all competitive programming fenwick tree performance tricks that there are, and put them
all in this crate.

# Naming

This library is called `ftree`, because the base data structure is `FenwickTree`.

# Changelog

See [CHANGELOG.md](https://github.com/brurucy/ftree/blob/master/CHANGELOG.md).
