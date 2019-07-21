#pragma once

#include <cstdint>
#include <vector>
#include <array>

template <typename T, int N>
struct decode_audio_buffer_t {
	std::size_t write_index;
	std::size_t read_index;
	std::array<T, N> buf;
    std::atomic<bool> paused;

	decode_audio_buffer_t() : 
		write_index(0),
		read_index(0),
        paused(false)
	{

	}

    void pause()
    {
        paused = true;
    }

    void resume()
    {
        paused = false;
    }

	bool write(const std::vector<T>& samples) 
	{
        if (paused)
        {
            return false;
        }
		else if ( (samples.size() + write_index) > buf.size() )
		{
			std::cout << "decode audio buffer full!" << std::endl;
			return false;
		}
		for (const auto& sample : samples) {
			buf[write_index] = sample;
			write_index++;
		}
		//std::cout << "Wrote " << samples.size() << " samples, " << buf.size() - write_index << " space remains" << std::endl;
		return true;
	}

	void clear() 
	{
		buf.fill(0.0f);
		write_index = 0;
		read_index = 0;
	}
};


