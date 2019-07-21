#pragma once

#include <atomic>
#include <cstdint>
#include <thread>
#include <chrono>

// Used internally as a single producer single consumer queue (ring buffer)
template <typename T>
struct ring_buffer_t {
	T* recs;
	// Use atomics for indices to prevent instruction reordering
	std::atomic_int read_index;
	std::atomic_int write_index;
	// Number of items the ring buffer can hold.
	size_t size;
	std::atomic<bool> terminateFlag;

	ring_buffer_t() :
		terminateFlag(false),
		recs(nullptr),
        size(0)
	{
	}

	bool initialize(const size_t sizeIn) 
    {
		read_index = 0;
		write_index = 0;
		size = sizeIn;
		recs = reinterpret_cast<T*>(malloc(sizeof(T) * size));
		return (recs != nullptr);
	}

	void terminate() 
    {
		terminateFlag = true;
	}

	~ring_buffer_t() 
	{
    /*
		if (nullptr != recs)
		{
			for (size_t k = 0; k < size; ++k) {
				free(recs[k]);
			}
			free(recs);
		}
        */
	}

	bool wait_for_empty_slot() const 
    {
		while ((read_index == write_index + 1) || (read_index == 0 && static_cast<int64_t>(write_index) == static_cast<int64_t>(size) - 1)) {
			std::this_thread::sleep_for(std::chrono::microseconds(100));
			if (terminateFlag) {
				return false;
			}
		}
		return true;
	}

	void inc_write_index()
    {
		// cast to signed so we don't break subtraction
		if (static_cast<int64_t>(write_index) == static_cast<int64_t>(size) - 1) {
			write_index = 0;
		}
		else {
			write_index++;
		}
	}

    T& pop_ref()
    {
        wait_for_data();
        T& curr = recs[read_index];
        if (read_index == size - 1) {
            read_index = 0;
        }
        else {
            read_index++;
        }
        return curr;
    }

	T pop() 
    {
		wait_for_data();
		T curr = recs[read_index];
		if (read_index == size - 1) {
			read_index = 0;
		}
		else {
			read_index++;
		}
		return curr;
	}

	bool wait_for_data() const 
    {
		while (empty()) {
			std::this_thread::sleep_for(std::chrono::microseconds(100));
			if (terminateFlag) {
				return false;
			}
		}
		return true;
	}

    bool empty() const {
        return read_index == write_index;
    }

};
