// Spreading independent jobs over the machine's processor cores.
#pragma once
#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>

namespace msv {

// Runs f(0) ... f(n-1), several at a time.
template <class F>
void parallel_for(int n, F f) {
    if (n <= 0) return;
    unsigned nt = std::min<unsigned>((unsigned)n, std::max(1u, std::thread::hardware_concurrency()));
    std::atomic<int> next(0);
    auto work = [&] {
        for (int k; (k = next++) < n;) f(k);
    };
    std::vector<std::thread> pool;
    for (unsigned i = 1; i < nt; i++) pool.emplace_back(work);
    work();
    for (auto& t : pool) t.join();
}

}  // namespace msv
