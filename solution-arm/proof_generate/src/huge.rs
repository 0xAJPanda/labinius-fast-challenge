
use std::alloc::{GlobalAlloc, Layout, System};

const HUGE_PAGE: usize = 2 << 20;
const MIN_ADVISED: usize = 2 * HUGE_PAGE;

pub struct HugePages;

#[cfg(target_os = "linux")]
fn advise(p: *mut u8, size: usize) {
    extern "C" {
        fn madvise(addr: *mut core::ffi::c_void, len: usize, advice: i32) -> i32;
    }
    const MADV_HUGEPAGE: i32 = 14;
    if p.is_null() || size < MIN_ADVISED {
        return;
    }
    let start = (p as usize).next_multiple_of(HUGE_PAGE);
    let end = (p as usize + size) & !(HUGE_PAGE - 1);
    if end > start {
        unsafe {
            madvise(start as *mut core::ffi::c_void, end - start, MADV_HUGEPAGE);
        }
    }
}

#[cfg(not(target_os = "linux"))]
fn advise(_p: *mut u8, _size: usize) {}

unsafe impl GlobalAlloc for HugePages {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        advise(p, layout.size());
        p
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc_zeroed(layout);
        advise(p, layout.size());
        p
    }

    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        System.dealloc(p, layout)
    }

    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let q = System.realloc(p, layout, size);
        advise(q, size);
        q
    }
}
