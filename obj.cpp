#include "obj.h"
#include "appx"

inline bool is_space(char c) { return c == ' ' || c == '\t'; }
inline bool is_endl(char c) { return c == '\r' || c == '\n'; }
inline bool is_sep(char c) { return is_space(c) || is_endl(c); }

Obj parse_object(char const* objText)
{
	Obj obj;

	size_t fileSize = strlen(objText);
	size_t percentDelta = (fileSize + 99) / 100;
	char const *nextPercentMarker = objText;

	appx::Task task("obj");

	bool normalsReordered = false;
	std::vector<glm::pod_vec<unsigned, 3>> nf;

	for (char const* objCursor = objText; *objCursor; )
	{
		// Trim
		while (is_space(*objCursor)) ++objCursor;

		switch (*objCursor)
		{
		case 'v':
			{
				++objCursor;

				switch (*objCursor)
				{
				// Position
				case '\t':
				case ' ':
					{
						++objCursor;
						
						glm::pod_vec<float, 3> vec;
						char* numEnd;
						vec.c[0] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						vec.c[1] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						vec.c[2] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						obj.v.push_back(vec);
					}
					break;
				// Normal
				case 'n':
					{
						++objCursor;

						glm::pod_vec<float, 3> vec;
						char* numEnd;
						vec.c[0] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						vec.c[1] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						vec.c[2] = (float) strtod(objCursor, &numEnd); objCursor = numEnd;
						obj.n.push_back(vec);
					}
					break;
				}
			}
			break;

		// Face
		case 'f':
			{
				++objCursor;

				glm::pod_vec<unsigned, 3> vind = { };
				glm::pod_vec<unsigned, 3> nind = { };

				auto nv = (long) obj.v.size();
				auto nn = (long) obj.n.size();

				for (int i = 0; i < 3; ++i)
				{
					char* numEnd;
					long l;
					
					l = strtol(objCursor, &numEnd, 10);
					if (l < 0) l += nv;
					else --l;
					vind.c[i] = (unsigned) l;
					objCursor = numEnd;

					if (*objCursor == '/')
					{
						++objCursor;
						l = strtol(objCursor, &numEnd, 10);
						objCursor = numEnd;

						if (*objCursor == '/')
						{
							++objCursor;
							l = strtol(objCursor, &numEnd, 10);
							if (l < 0) l += nn;
							else --l;
							nind.c[i] = (unsigned) l;
							objCursor = numEnd;

							normalsReordered |= nind.c[i] != vind.c[i];
						}
					}
				}

				obj.f.push_back(vind);
				nf.push_back(nind);
			}
			break;

		case 'u':
			{
				if (strncmp(objCursor, "usemtl", arraylen("usemtl") - 1) == 0)
				{
					objCursor += arraylen("usemtl") - 1;

					// Extract name
					while (is_space(*objCursor)) ++objCursor;
					auto mtlNameStart = objCursor;
					while (*objCursor && !is_endl(*objCursor)) ++objCursor;
					auto mtlNameEnd = objCursor;

					obj.m.emplace_back((unsigned) obj.f.size(), std::string(mtlNameStart, mtlNameEnd));
				}
			}
			break;
		}

		// Next line
		while (*objCursor && !is_endl(*objCursor)) ++objCursor;
		while (is_endl(*objCursor)) ++objCursor;

		if (objCursor >= nextPercentMarker)
		{
			auto pct = 100 * (objCursor - objText) / fileSize;
			task.progressPct((float) pct);
			nextPercentMarker = objText + std::min((pct + 10) * percentDelta, fileSize);
		}
	}

	if (normalsReordered)
	{
		auto unorderedN = std::move(obj.n);
		obj.n.resize(obj.v.size());

		for (size_t i = 0, ie = nf.size(); i < ie; ++i)
			for (size_t j = 0; j < 3; ++j)
				obj.n[obj.f[i].c[j]] = unorderedN[nf[i].c[j]];
	}

	return obj;
}
