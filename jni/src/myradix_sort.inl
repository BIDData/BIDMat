/*
 *  Copyright 2008-2012 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <thrust/detail/config.h>

// do not attempt to compile this file with any other compiler
#if THRUST_DEVICE_COMPILER == THRUST_DEVICE_COMPILER_NVCC

#include <thrust/detail/copy.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/iterator/iterator_traits.h>

#include <thrust/detail/uninitialized_array.h>
#include <thrust/detail/type_traits.h>
#include <thrust/detail/util/align.h>


__THRUST_DISABLE_MSVC_POSSIBLE_LOSS_OF_DATA_WARNING_BEGIN


#include <thrust/detail/backend/cuda/detail/b40c/radixsort_api.h>

///////////////////////
// Key-Value Sorting //
///////////////////////

template<typename RandomAccessIterator1,
         typename RandomAccessIterator2>
int radix_sort_spine(RandomAccessIterator1 first1,
                     RandomAccessIterator1 last1,
                     RandomAccessIterator2 first2)
{
    typedef typename thrust::iterator_value<RandomAccessIterator1>::type K;
    typedef typename thrust::iterator_value<RandomAccessIterator2>::type V;
    unsigned int num_elements = last1 - first1;
    thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<K,V> sorter(num_elements);
    return sorter.SpineElements();
}

// sort values directly
template<typename RandomAccessIterator1,
         typename RandomAccessIterator2>
void stable_radix_sort_by_key(RandomAccessIterator1 first1,
                              RandomAccessIterator1 last1,
                              RandomAccessIterator2 first2,
                              RandomAccessIterator1 first3,
                              RandomAccessIterator2 first4,
                              int * d_spine,
                              bool * d_alt_storage)
{
  //    typedef thrust::detail::cuda_device_space_tag space;
    typedef typename thrust::iterator_value<RandomAccessIterator1>::type K;
    typedef typename thrust::iterator_value<RandomAccessIterator2>::type V;
    
    unsigned int num_elements = last1 - first1;

    thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<K,V> sorter(num_elements);
    thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortStorage<K,V>    storage;

    // define storage
    storage.d_keys             = thrust::raw_pointer_cast(&*first1);
    storage.d_values           = thrust::raw_pointer_cast(&*first2);
    storage.d_alt_keys         = thrust::raw_pointer_cast(&*first3);
    storage.d_alt_values       = thrust::raw_pointer_cast(&*first4);
    storage.d_spine            = d_spine;
    storage.d_from_alt_storage = d_alt_storage;

    // perform the sort
    sorter.EnactSort(storage);
    
    // radix sort sometimes leaves results in the alternate buffers
    if (storage.using_alternate_storage)
    {
    printf("alternative storage\n");
        cudaMemcpy(storage.d_keys, storage.d_alt_keys, num_elements*sizeof(K), cudaMemcpyDeviceToDevice);
        cudaMemcpy(storage.d_values, storage.d_alt_values, num_elements*sizeof(V), cudaMemcpyDeviceToDevice);
    }
}


__THRUST_DISABLE_MSVC_POSSIBLE_LOSS_OF_DATA_WARNING_END


#endif // THRUST_DEVICE_COMPILER == THRUST_DEVICE_COMPILER_NVCC

