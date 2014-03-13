#pragma once

#include "stdx"
#include "mathx"
#include <vector>

struct Obj
{
	MOVE_GENERATE(Obj, MOVE_5, MEMBER, v, MEMBER, n, MEMBER, t, MEMBER, f, MEMBER, m);

	std::vector<glm::pod_vec<float, 3>> v;
	std::vector<glm::pod_vec<float, 3>> n;
	std::vector<glm::pod_vec<float, 2>> t;
	std::vector<glm::pod_vec<unsigned, 3>> f;
	std::vector<std::pair<unsigned, std::string>> m;

	Obj() { }
};

Obj parse_object(char const* objText);