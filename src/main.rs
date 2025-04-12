fn main() {
    for i in [3, 4, 5, 6, 8, 10, 12, 14] {
        let start = std::time::Instant::now();
        let ds = ff_sbibd::find_diffset(i - 1);
        println!(
            "Diffset order {} found in {:?} = {:?}",
            i - 1,
            start.elapsed(),
            ds
        );
    }
}
