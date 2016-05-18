
#ifndef WANG_RAND

uint tea(uint val0, uint val1, uint backoff = 16)
{
	uint v0 = val0;
	uint v1 = val1;
	uint s0 = 0;

	for (uint n = 0; n < backoff; n++)
	{
		s0 += 0x9e3779b9;
		v0 += ((v1<<4)+0xa341316c)^(v1+s0)^((v1>>5)+0xc8013ea4);
		v1 += ((v0<<4)+0xad90777d)^(v0+s0)^((v0>>5)+0x7e95761e);
	}

	return v0;
}

uint lcg(inout uint prev)
{
	const uint LCG_A = 1664525u;
	const uint LCG_C = 1013904223u;
	prev = (LCG_A * prev + LCG_C);
	return prev & 0x00FFFFFF;
}

uint lcg2(inout uint prev)
{
	prev = (prev*8121 + 28411)  % 134456;
	return prev;
}

// generate random float in [0, 1)
float rnd(inout uint prev)
{
	return float(lcg(prev)) / float(0x01000000);
}

// generate -1.0f or 1.0f
float rnd_sgn(inout uint prev)
{
	return (lcg(prev) >= 0x00800000) ? 1.0f : -1.0f;
}

struct RandState
{
	uint prev;
};

RandState rand_init(uint val0, uint val1, uint stride_unused = 0, uint backoff = 16)
{
	RandState state = { tea(val0, val1, backoff) };
	return state;
}

RandState rand_renew(RandState oldState, uint val1, uint backoff = 16)
{
	RandState state = oldState;
	state.prev = tea(oldState.prev, val1, backoff);
	return state;
}

float rand_next(inout RandState s)
{
	return rnd(s.prev);
}

float rand_next_sign(inout RandState s)
{
	return rnd_sgn(s.prev);
}

#else

// Random numbers: http://www.reedbeta.com/blog/2013/01/12/quick-and-easy-gpu-random-numbers-in-d3d11/

uint wang_hash(uint seed)
{
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
}

uint wang_init(uint val0, uint seed)
{
	return 3U * val0 + seed * 73495523U;
}

uint wang_sample_hash(inout uint state, uint stride)
{
	uint hash = wang_hash(state);
	state += stride;
	return hash;
}
float wang_sample(inout uint state, uint stride)
{
//	return float(wang_sample_hash(state) & 0xFFFFF) / float(0xFFFFF);
	return wang_sample_hash(state, stride) / 4294967296.0f;
}
float wang_sample_sign(inout uint state, uint stride)
{
	return (wang_sample_hash(state, stride) <= 0x7fffffffU) ? 1.0f : -1.0f;
}

struct RandState
{
	uint state;
	uint stride;
};

RandState rand_init(uint val0, uint val1, uint stride, uint backoff_unused = 0)
{
	RandState state;
	state.state = wang_init(val0, val1);
	state.stride = stride;
	return state;
}

RandState rand_renew(RandState oldState, uint val1, uint backoff = 16)
{
	RandState state = oldState;
	state.prev = wang_init(oldState.prev, val1);
	return state;
}

float rand_next(inout RandState s)
{
	return wang_sample(s.state, s.stride);
}

float rand_next_sign(inout RandState s)
{
	return wang_sample_sign(s.state, s.stride);
}


#endif
