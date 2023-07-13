#![cfg_attr(not(test), no_std)]
extern crate alloc;
use alloc::vec::Vec;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

const fn least_significant_bit(idx: usize) -> usize {
    return idx & idx.wrapping_neg();
}

const fn most_significant_bit(idx: usize) -> usize {
    if idx == 0 {
        return 0;
    }
    return 1 << (usize::BITS - 1 - idx.leading_zeros());
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, PartialEq)]
pub struct FenwickTree {
    inner: Vec<usize>,
}

impl FromIterator<usize> for FenwickTree {
    fn from_iter<T: IntoIterator<Item=usize>>(iter: T) -> Self {
        return FenwickTree::new(iter, |val| *val)
    }
}

impl<const N: usize> From<[usize; N]> for FenwickTree
{
    fn from(value: [usize; N]) -> Self {
        return FenwickTree::from_iter(value.into_iter());
    }
}

impl FenwickTree {
    /// Creates a new fenwick tree.
    ///
    /// # Examples
    ///
    /// ```
    /// use ftree::FenwickTree;
    ///
    /// let lengths: [usize; 5] = [1, 6, 3, 9, 2];
    /// // This is how lengths fenwick tree will look like internally
    /// let fenwick_tree: Vec<usize> = vec![1, 7, 3, 19, 2];
    /// // And this is how it can be constructed with new
    /// let initialized_fenwick_tree: FenwickTree = FenwickTree::new(lengths, |length| *length);
    /// ```
    pub fn new<T>(container: impl IntoIterator<Item = T>, f: impl Fn(&T) -> usize) -> Self {
        let mut inner: Vec<usize> = container.into_iter().map(|value| f(&value)).collect();
        let n = inner.len();

        for i in 0..n {
            let parent = i | (i + 1);
            if parent < n {
                inner[parent] += inner[i];
            }
        }

        FenwickTree { inner }
    }
    /// Computes the prefix sum up until the desired index.
    ///
    /// The prefix sum up until the zeroth element is 0, since there is nothing before it.
    ///
    /// The prefix sum up until an index larger than the length is undefined, since every
    /// element after the length - 1 is undefined.
    ///
    /// # Examples
    ///
    /// ```
    /// use ftree::FenwickTree;
    ///
    /// let lengths = [1, 6, 3, 9, 2];
    /// let fenwick_array = FenwickTree::new(lengths, |item| *item);
    ///
    /// let cases: Vec<(usize, usize)> =
    ///  vec![(0, 0), (1, 1), (2, 7), (3, 10), (4, 19), (5, 21), (6, 0)];
    ///
    /// cases
    ///   .into_iter()
    ///   .for_each(|(idx, expected_sum)| assert_eq!(fenwick_array.prefix_sum(idx), expected_sum))
    /// ```
    pub fn prefix_sum(&self, index: usize) -> usize {
        let mut sum = 0;
        let mut current_idx = index;

        if index > self.inner.len() {
            return sum;
        }

        while current_idx > 0 {
            sum += self.inner[current_idx - 1];
            current_idx &= current_idx - 1
        }

        return sum;
    }
    /// Either increments, or decrements, the fenwick tree at a given index. Take care to not render
    /// it invalid.
    ///
    /// # Examples
    ///
    /// ```
    /// use ftree::FenwickTree;
    ///
    /// let lengths = [1, 6, 3, 9, 2];
    /// let mut fenwick_array = FenwickTree::new(lengths, |item| *item);
    ///
    /// let cases: Vec<(usize, usize)> = vec![(0, 0), (1, 2), (2, 8), (3, 11), (4, 20), (5, 22), (6, 0)];
    ///
    /// fenwick_array.update_index(0, 1);
    ///
    /// cases
    ///   .into_iter()
    ///   .for_each(|(idx, expected_sum)| assert_eq!(fenwick_array.prefix_sum(idx), expected_sum))
    /// ```
    pub fn update_index(&mut self, index: usize, diff: isize) {
        let mut current_idx = index;
        while current_idx < self.inner.len() {
            self.inner[current_idx] = (self.inner[current_idx] as isize + diff) as usize;
            current_idx |= current_idx + 1
        }
    }
    /// Given a sum, finds the slot in which in which it would be "contained" within the original
    /// array.
    ///
    /// # Examples
    ///
    /// ```
    /// use ftree::FenwickTree;
    ///
    /// let lengths = [1, 6, 3, 9, 2];
    /// let mut fenwick_array = FenwickTree::new(lengths, |item| *item);
    ///
    /// let cases: Vec<(usize, usize)> = vec![(0, 0), (6, 1), (9, 2), (18, 3), (20, 4)];
    ///
    /// cases
    ///   .into_iter()
    ///   .for_each(|(prefix_sum, idx)| assert_eq!(fenwick_array.index_of(prefix_sum), idx))
    /// ```
    pub fn index_of(&self, prefix_sum: usize) -> usize {
        let length = self.inner.len();
        let mut prefix_sum = prefix_sum;
        let mut idx = 0;
        let mut x = most_significant_bit(length) * 2;
        while x > 0 {
            let lsb = least_significant_bit(x);
            let half_lsb = lsb / 2;
            if x <= length && self.inner[x - 1] < prefix_sum {
                idx = x;
                prefix_sum -= self.inner[x - 1];
                x += half_lsb;
            } else {
                if lsb % 2 > 0 {
                    break;
                }
                x = x + half_lsb - lsb;
            }
        }
        idx
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use alloc::vec::Vec;
    use super::FenwickTree;

    #[test]
    fn test_new() {
        let lengths: [usize; 5] = [1, 6, 3, 9, 2];
        let expected_index: Vec<usize> = vec![1, 7, 3, 19, 2];
        let actual_index = FenwickTree::new(lengths, |item| *item);
        assert_eq!(expected_index, actual_index.inner)
    }

    #[test]
    fn test_prefix_sum() {
        let lengths = [1, 6, 3, 9, 2];
        let fenwick_array = FenwickTree::new(lengths, |item| *item);

        let cases: Vec<(usize, usize)> =
            vec![(0, 0), (1, 1), (2, 7), (3, 10), (4, 19), (5, 21), (6, 0)];
        // The prefix sum up until the zeroth element is 0, since there is nothing before it
        // The prefix sum up until an index larger than the length is undefined, since every
        // element after the length - 1 is undefined
        cases
            .into_iter()
            .for_each(|(idx, expected_sum)| assert_eq!(fenwick_array.prefix_sum(idx), expected_sum))
    }

    #[test]
    fn test_update_index() {
        let lengths = [1, 6, 3, 9, 2];
        let mut fenwick_array = FenwickTree::new(lengths, |item| *item);

        let cases: Vec<(usize, usize)> = vec![(0, 2), (1, 8), (2, 3), (3, 20), (4, 2)];

        fenwick_array.update_index(0, 1);

        cases
            .into_iter()
            .for_each(|(idx, expected_value)| assert_eq!(fenwick_array.inner[idx], expected_value))
    }

    #[test]
    fn test_index_of() {
        let lengths: Vec<usize> = vec![1, 6, 3, 9, 2];
        let fenwick_array = FenwickTree::new(lengths, |item| *item);

        let cases: Vec<(usize, usize)> = vec![(0, 0), (6, 1), (9, 2), (18, 3), (20, 4)];

        cases
            .into_iter()
            .for_each(|(prefix_sum, idx)| assert_eq!(fenwick_array.index_of(prefix_sum), idx))
    }
}
