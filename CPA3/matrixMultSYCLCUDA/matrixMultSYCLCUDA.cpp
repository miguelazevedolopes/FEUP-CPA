#include <CL/sycl.hpp>

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>

using namespace cl::sycl;

class mxm_kernel;

template <typename T>
bool local_mxm(cl::sycl::queue &q, T *MA, T *MB, T *MC, int matSize, int blockSize)
{
  {
    /* Buffers can be constructed with property lists. In this example,
     * the buffer is given the property "use host pointer", which tells
     * the runtime to use the host pointer for all data storage (instead
     * of making copies internally). Additionally, when running on a
     * device that shares memory with the host (for example a CPU),
     * "zero-copy" memory optimisations can be used by the driver. */
    range<1> dimensions(matSize * matSize);
    const property_list props = {property::buffer::use_host_ptr()};
    buffer<T> bA(MA, dimensions, props);
    buffer<T> bB(MB, dimensions, props);
    buffer<T> bC(MC, dimensions, props);

    q.submit([&](handler &cgh)
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
  }
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
      MC[i * matSize + j] = 0.0f; // i * matSize + j;
    }
  }
}

int main(int argc, char *argv[])
{
  float *MA;
  float *MB;
  float *MC;

  for (int n = 1024; n <= 8192; n += 1024) {
    MA = new float[n * n];
    MB = new float[n * n];
    MC = new float[n * n];

    std::cout << "Matrix Size N = " << n << std::endl
              << std::endl;

    queue q(gpu_selector_v);

    initMatrix(MA, MB, MC, n);
    std::cout << "CUDA with " << q.get_device().get_info<sycl::info::device::name>() << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    local_mxm(q, MA, MB, MC, n, 32);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl
              << std::endl;

    delete[] MA;
    delete[] MB;
    delete[] MC;
  }

  return 0;
}
