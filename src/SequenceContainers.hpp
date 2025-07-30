//  SequenceContainers.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3

/**
 * DynamicBufferArray
 *
 * DynamicBufferArray is a templated dynamic array class that
 * provides similar functionalities to the C++ Standard Library's std::vector.
 * It manages its memory internally and resizes the array dynamically, offering
 * better memory management compared to traditional C++ arrays.
 *
 * The class uses a fixed-size buffer for the first `fixed_capacity` elements
 * to avoid heap allocations until necessary. The class adheres to RAII principles.
 */


#pragma once
#include "ofMain.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>




template <typename T>
class DynamicBufferArray
{
public:
	
	// ------------- Constructors and Destructor -------------
	DynamicBufferArray();
	DynamicBufferArray(size_t _capacity) ;
	DynamicBufferArray(const DynamicBufferArray& other);
	DynamicBufferArray& operator=(const DynamicBufferArray& other);
	~DynamicBufferArray();
	
	
	// ------------- Element Access -------------
	T& operator[](size_t index); //must return as reference as array element can be put on left side
	const T& operator[](size_t index) const;
	
	
	// ------------- Capacity and Size -------------
	bool empty() const; // Returns size() == 0.
	int size() const;  // Returns the number of elements in the list.
	
	
	void reserve(size_t _size); //Reserves a minimum capacity for the array.
	void resize(size_t _size, const T& fill); //resizes the list to contain _size elements, filling new elements with fill
	
	
	
	// ------------- Modifiers -------------
	void push_back(const T& _data);
	T pop_back();
	void clear();
	void swap(DynamicBufferArray &other);
	
	template <typename... Args>
	void emplace_back(Args&&... args);
	
	
private:
	// ------------- Internal storage and metadata for DynamicBufferArray. -------------
	
	enum { fixed_cap = 64 };
	struct alignas(16) ListData //is memory-aligned for optimal access speed. In C++11 and newer, you can use the alignas keyword.
	{
	ListData();
	T buffer[fixed_cap];// Stores a fixed-size buffer in advance to avoid requiring a heap allocation until we run out of space.
	T* data;  // Pointer to the actual data
	size_t size;   //total current size of array/vector, i.e., number of elements
	size_t capacity;  // Total capacity of the allocated array, i.e., max size
};
ListData listData;  // Internal data storage


};



template <typename T>
DynamicBufferArray<T>::ListData::ListData() : data(buffer), size(0), capacity(fixed_cap){}


template <typename T>
DynamicBufferArray<T>::DynamicBufferArray() {}


template <typename T>
DynamicBufferArray<T>::DynamicBufferArray(const DynamicBufferArray& other)
	{
	if (other.listData.capacity == fixed_cap)
	{
		listData = other.listData;
		listData.data = listData.buffer;
	}
	else
	{
		reserve(other.listData.size);
		for (int j=0; j < other.size(); ++j)
		{
			listData.data[j] = other.listData.data[j];
		}
		listData.size = other.listData.size;
		listData.capacity = other.listData.capacity;
	}
}




template <typename T>
DynamicBufferArray<T> &DynamicBufferArray<T>::operator=(const DynamicBufferArray<T>& other)
	{
	DynamicBufferArray(other).swap(*this);
	return *this;
}




template <typename T>
DynamicBufferArray<T>::~DynamicBufferArray()
	{
	if(listData.data != listData.buffer)
	{
		delete[] listData.data;
	}
}



// Returns size() == 0.
template <typename T>
bool DynamicBufferArray<T> :: empty() const
	{
	return listData.size == 0;
}

// Returns the number of elements in the list.
template <typename T>
int DynamicBufferArray<T> :: size() const
	{
	return listData.size;
}

template <typename T>
T &DynamicBufferArray<T>::operator[](size_t index)
	{
	assert(index >= 0 && index < listData.size);
	return listData.data[index];
}

template <typename T>
const T &DynamicBufferArray<T>::operator[](size_t index) const
	{
	assert(index >= 0 && index < listData.size);
	return listData.data[index];
}


template <typename T>
void DynamicBufferArray<T> ::clear()
	{
	listData.size = 0;
}





template <typename T>
void DynamicBufferArray<T> ::reserve(size_t _size)
{
	if (_size > listData.capacity)
	{
		T* newMemory = new T[_size];
		if (!newMemory)
		{
			throw std::bad_alloc();
		}
		std::copy(listData.data, listData.data + listData.capacity, newMemory);
		if (listData.data != listData.buffer)
		{
			delete[] listData.data;
		}
		
		listData.data = newMemory;
		listData.capacity = _size;
	}
}




template <typename T>
void DynamicBufferArray<T> ::resize(size_t _size, const T &fill)
{
	if (_size > listData.capacity)
	{
		reserve(_size + 1);
	}
	if (_size > listData.size)
	{
		for (int j = listData.size; j < _size; ++j)
		{
			listData.data[j] = fill;
		}
	}
	listData.size = _size;
}








template <typename T>
void DynamicBufferArray<T> :: push_back(const T &_data)
{
	if (listData.size >= listData.capacity)
	{
		//reserve(listData.capacity * 2);
		reserve(static_cast<size_t>(listData.capacity * 1.618)); // Golden Ratio
	}
	listData.data[listData.size++] = _data;
	
}



template <typename T>
T DynamicBufferArray<T> :: pop_back()
{
	return listData.data[--listData.size];
	
}








template <typename T>
void DynamicBufferArray<T>::swap(DynamicBufferArray &other)
{
	ListData& ld1 = listData;
	ListData& ld2 = other.listData;
	
	const int use_fixed1 = ld1.data == ld1.buffer;
	const int use_fixed2 = ld2.data == ld2.buffer;
	
	const ListData temp = ld1;
	ld1 = ld2;
	ld2 = temp;
	
	if (use_fixed1)
		ld2.data = ld2.buffer;
	if (use_fixed2)
		ld1.data = ld1.buffer;
}






template <typename T>
template <typename... Args>
void DynamicBufferArray<T>::emplace_back(Args&&... args)
{
	// Make sure there's enough space for a new element
	if (listData.size >= listData.capacity)
	{
		reserve(listData.capacity * 2);
	}
	// Construct the object in-place at the end of the array using std::forward to preserve argument types
	new (&listData.data[listData.size]) T(std::forward<Args>(args)...);
	++listData.size;
}



