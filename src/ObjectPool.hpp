//  ObjectPool.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3


#pragma once
#include "ofMain.h"
#include <atomic>
#include <future>
#include <thread>
#include <functional>
#include <stdexcept>




template <class T>
class ObjectPool
{
public:
	ObjectPool() {}
	
	ObjectPool(size_t size)
	{
		for (size_t i = 0; i < size; ++i)
		{
			auto object = new T();
			pool.push(object);
		}
	}
	
	void releaseAll()
	{
		while (!pool.empty())
		{
			delete pool.top();
			pool.pop();
		}
	}
	
	
	~ObjectPool()
	{
		releaseAll();
	}
	
	T* acquire()
	{
		if (pool.empty())
		{
			// Consider either dynamically growing the pool or return nullptr.
			// Here's an example of growing the pool:
			auto object = new T();
			pool.push(object);
			
			//std::cout << "\n\nThe object pool is empty";
			//return nullptr; // if you don't want to grow the pool
		}
		
		auto obj = pool.top();
		pool.pop();
		return obj;
	}
	
	
	
	void adjustPoolSize(size_t newSize)
	{
		size_t currentSize = pool.size();
		
		// If the new size is greater than the current size, allocate more objects
		if (newSize > currentSize)
		{
			allocate(newSize - currentSize);
		}
		// If the new size is smaller than the current size, delete extra objects
		else if (newSize < currentSize)
		{
			for (size_t i = 0; i < (currentSize - newSize); ++i)
			{
				delete pool.top();
				pool.pop();
			}
		}
		
	}
	void release(T* obj) //reset the state of the objects before you reinsert them into the pool.
	{
		if (obj) // ensure the object is not nullptr
		{
			obj->reset();
			pool.push(obj);
		}
	}
	
	size_t size()
	{
		return(pool.size());
	}
	
	
private:
	std::stack<T*> pool; //data structure of available/reusable pre-allocated  objects. efficient way to add objects (when they are created or returned), and to remove objects (when they are acquired) in constant time, O(1).
	
	void allocate(size_t _size)
	{
		for (size_t i = 0; i < _size; ++i)
		{
			auto object = new T();
			pool.push(object);
		}
		
	}
	
	
};




