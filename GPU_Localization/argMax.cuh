#include <cub/cub.cuh>
#define WARP_SIZE 32
template <typename K, typename V, unsigned int TPB>
  class BlockArgMax {

    using Reducer =
      cub::BlockReduce<cub::KeyValuePair<K, V>, TPB>;

  public:

    //
    // @brief compute argmax, the key associated with the largest value
    // @param key key element of each thread's pair
    // @param value value element of each thread's pair
    // @return key associated with largest value (in 0th thread ONLY)
    //
    __device__
    static K argmax(const K &key, const V &value)
    {
      cub::KeyValuePair<K, V> myPair;
      myPair.key   = key;
      myPair.value = value;

      __shared__ typename Reducer::TempStorage CUB_tmp;

      auto resultPair =
        Reducer(CUB_tmp).Reduce(myPair, cub::ArgMax());

      return resultPair.key;
    }

    //
    // @brief compute argmax, the key associated with the largest value,
 //    in a way that supports ops over arrays greater than block size
    //
    // @param key key element of each thread's pair
    // @param value value element of each thread's pair
    // @param maxValue output parameter for largest value (set in 0th
    //          thread ONLY)
    // @param nThreads number of threads with valid inputs
    // @return key associated with largest value (in 0th thread ONLY)
    //
    __device__
    static K argmax(const K &key, const V &value, V &maxValue,
                    unsigned int nThreads = TPB)
    {
      cub::KeyValuePair<K, V> myPair;
      myPair.key   = key;
      myPair.value = value;

      __shared__ typename Reducer::TempStorage CUB_tmp;

      auto resultPair =
        Reducer(CUB_tmp).Reduce(myPair, cub::ArgMax(), nThreads);

      if (threadIdx.x == 0)
        maxValue = resultPair.value;

      return resultPair.key;
    }

  };

  template <typename T, unsigned int TPB>
  __device__
  T broadcast(const T &v, unsigned int sourceTID = 0)
  {
    if (TPB == WARP_SIZE)
      {
        return __shfl_sync((1ULL << WARP_SIZE) - 1, v, sourceTID);
      }
    else
      {
        __shared__ T vShared;

        __syncthreads();
        if (threadIdx.x == sourceTID)
          vShared = v;
        __syncthreads();

        return vShared;
      }
  }

  template <typename T, unsigned int TPB>
  class BlockReduce {

    using Reducer =
      cub::BlockReduce<T, TPB>;

  public:

    //
    // @brief compute a sum reduction over the first nThreads threads
    //   of the block
    // @param v values to sum
    // @param nThreads # of threads to sum over (defaults to all)
    //
    // @returns sum (in thread 0 ONLY)
    //
    __device__
    static T sum(const T &v, unsigned int nThreads = TPB)
    {
      __shared__ typename Reducer::TempStorage CUB_tmp;

      return Reducer(CUB_tmp).Sum(v, nThreads);
    }

    __device__
    static T max(const T &v, unsigned int nThreads = TPB)
    {
      __shared__ typename Reducer::TempStorage CUB_tmp;

      return Reducer(CUB_tmp).Reduce(v, cub::Max());
    }

    __device__
    static T min(const T &v, unsigned int nThreads = TPB)
    {
      __shared__ typename Reducer::TempStorage CUB_tmp;

      return Reducer(CUB_tmp).Reduce(v, cub::Min());
    }
  };
