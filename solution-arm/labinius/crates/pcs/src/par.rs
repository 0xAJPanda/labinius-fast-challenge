
use std::sync::OnceLock;

pub fn threads() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| {
        std::env::var("LABINIUS_THREADS")
            .ok()
            .and_then(|s| s.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
    })
}

#[inline]
pub fn block(n: usize, parts: usize, t: usize, align: usize) -> (usize, usize) {
    let units = n / align;
    let lo = units * t / parts;
    let hi = units * (t + 1) / parts;
    let end = if t + 1 == parts { n } else { hi * align };
    (lo * align, end)
}

pub fn split_mut<'a, T>(s: &'a mut [T], edges: &[(usize, usize)]) -> Vec<&'a mut [T]> {
    let mut out = Vec::with_capacity(edges.len());
    let mut rest = s;
    let mut base = 0;
    for &(lo, hi) in edges {
        let (_, tail) = core::mem::take(&mut rest).split_at_mut(lo - base);
        let (mid, tail) = tail.split_at_mut(hi - lo);
        out.push(mid);
        rest = tail;
        base = hi;
    }
    out
}

#[inline]
pub fn workers(n: usize, min: usize) -> usize {
    (n / min.max(1)).clamp(1, threads())
}

pub fn map_blocks<R: Send>(
    n: usize,
    parts: usize,
    align: usize,
    f: impl Fn(usize, usize, usize) -> R + Sync,
) -> Vec<R> {
    if parts <= 1 {
        return vec![f(0, 0, n)];
    }
    let f = &f;
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..parts - 1)
            .map(|t| {
                let (s, e) = block(n, parts, t, align);
                scope.spawn(move || f(t, s, e))
            })
            .collect();
        let (s, e) = block(n, parts, parts - 1, align);
        let last = f(parts - 1, s, e);
        let mut out: Vec<R> = handles
            .into_iter()
            .map(|h| h.join().expect("worker panicked"))
            .collect();
        out.push(last);
        out
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blocks_tile_the_range() {
        for n in [0usize, 1, 7, 8, 63, 64, 65, 1000, 4096] {
            for parts in 1..=9 {
                for align in [1usize, 8, 64] {
                    let mut next = 0;
                    for t in 0..parts {
                        let (s, e) = block(n, parts, t, align);
                        assert_eq!(s, next);
                        assert!(e >= s);
                        if t + 1 < parts {
                            assert_eq!(e % align, 0);
                        }
                        next = e;
                    }
                    assert_eq!(next, n);
                }
            }
        }
    }
}
