#include <CL/sycl.hpp>

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace cl::sycl;
#define NSEC_IN_MSEC 1000000.0

class mxm_kernel;

template <typename T>
bool local_mxm(cl::sycl::queue &q, T *MA, T *MB, T *MC, int matSize, int blockSize)
{
  range<1> dimensions(matSize * matSize);
  const property_list props = {property::buffer::use_host_ptr()};
  buffer<T> bA(MA, dimensions, props);
  buffer<T> bB(MB, dimensions, props);
  buffer<T> bC(MC, dimensions, props);

  sycl::event event = q.submit([&](handler &cgh)
           {
            auto pA = bA.template get_access<access::mode::read>(cgh);
            auto pB = bB.template get_access<access::mode::read>(cgh);
            auto pC = bC.template get_access<access::mode::write>(cgh);
            auto localRange = range<1>(blockSize * blockSize);

            accessor<T, 1, access::mode::read_write, access::target::local> pBA(
                localRange, cgh);
            accessor<T, 1, access::mode::read_write, access::target::local> pBB(
                localRange, cgh);

            cgh.parallel_for<mxm_kernel>(
                nd_range<2>{range<2>(matSize, matSize),
                range<2>(blockSize, blockSize)},
                [=](nd_item<2> it) {
                    // Current block
                    int blockX = it.get_group(1);
                    int blockY = it.get_group(0);

                    // Current local item
                    int localX = it.get_local_id(1);
                    int localY = it.get_local_id(0);

                    // Start in the A matrix
                    int a_start = matSize * blockSize * blockY;
                    // End in the b matrix
                    int a_end = a_start + matSize - 1;
                    // Start in the b matrix
                    int b_start = blockSize * blockX;

                    // Result for the current C(i,j) element
                    T tmp = 0.0f;
                    // We go through all a, b blocks
                    for (int a = a_start, b = b_start; a <= a_end;
                        a += blockSize, b += (blockSize * matSize)) {
                        // Copy the values in shared memory collectively
                        pBA[localY * blockSize + localX] =
                            pA[a + matSize * localY + localX];
                        // Note the swap of X/Y to maintain contiguous access
                        pBB[localX * blockSize + localY] =
                            pB[b + matSize * localY + localX];
                        it.barrier(access::fence_space::local_space);
                        // Now each thread adds the value of its sum
                        for (int k = 0; k < blockSize; k++) {
                            tmp +=
                                pBA[localY * blockSize + k] * pBB[localX * blockSize + k];
                        }
                        // The barrier ensures that all threads have written to local
                        // memory before continuing
                        it.barrier(access::fence_space::local_space);
                    }
                    auto elemIndex = it.get_global_id(0) * it.get_global_range()[1] +
                        it.get_global_id(1);
                    // Each thread updates its position
                    pC[elemIndex] = tmp;
                }); });

  event.wait();
  uint64_t start =
      event.get_profiling_info<sycl::info::event_profiling::command_start>();
  uint64_t end =
      event.get_profiling_info<sycl::info::event_profiling::command_end>();
  double duration = static_cast<double>(end - start) / NSEC_IN_MSEC;
  std::cout << "Time = " << std::fixed << std::setprecision(3) << duration << " msec" << std::endl << std::endl;

  return false;
}

void initMatrix(float *MA, float *MB, float *MC, int matSize)
{
  // Matrix initialization
  for (int i = 0; i < matSize; i++)
  {
    for (int j = 0; j < matSize; j++)
    {
      MA[i * matSize + j] = 0.0f;
      if (i == j)
      {
        MA[i * matSize + j] = 1.0f;
      }
      MB[i * matSize + j] = 2.0f;
      MC[i * matSize + j] = 0.0f;
    }
  }
}

int main(int argc, char *argv[])
{
  float *MA;
  float *MB;
  float *MC;

  sycl::property_list prop_list{sycl::property::queue::enable_profiling()};

  for (int n = 1024; n <= 8192; n += 1024)
  {
    MA = new float[n * n];
    MB = new float[n * n];
    MC = new float[n * n];

    std::cout << "Matrix Size N = " << n << std::endl;

    queue q(gpu_selector_v, prop_list);

    initMatrix(MA, MB, MC, n);
    std::cout << "CUDA with " << q.get_device().get_info<sycl::info::device::name>() << std::endl;
    local_mxm(q, MA, MB, MC, n, 32);

    delete[] MA;
    delete[] MB;
    delete[] MC;
  }

  return 0;
}
