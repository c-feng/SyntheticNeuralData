#include "mex.h"

#include <unistd.h>
#include <complex>
#include <map>
#include "mhs_error.h"
#include "mhs_vector.h"
#include <algorithm>

#define _SUPPORT_MATLAB_
#include "sta_mex_helpfunc.h"
// #define _SUPPORT_MATLAB_

#include <complex>
#include <vector>
#include <string>
#include <ctime>
#include <list>
#include <sstream>
#include <string>
#include <limits>
#include <omp.h>

#define SUB2IND(X, Y, Z, shape) (((Z) * (shape[1]) + (Y)) * (shape[2]) + (X))
#include "mhs_vector.h"

//mex mhs_distmapC.cc -lgomp CXXFLAGS=" -O3   -Wfatal-errors  -std=c++11 -fopenmp-simd -fopenmp  -fPIC -march=native"

template <typename T>
void hsv2rgb(T h, T s, T v, T &r, T &g, T &b)
{
	T hh, p, q, t, ff;
	int i;

	if (s <= 0.0)
	{
		r = v;
		g = v;
		b = v;
	}
	hh = h;
	if (hh >= 360.0)
	{
		hh = 0.0;
	}
	hh /= 60.0;
	i = (int)hh;
	ff = hh - i;
	p = v * (1.0 - s);
	q = v * (1.0 - (s * ff));
	t = v * (1.0 - (s * (1.0 - ff)));

	switch (i)
	{
	case 0:
		r = v;
		g = t;
		b = p;
		break;
	case 1:
		r = q;
		g = v;
		b = p;
		break;
	case 2:
		r = p;
		g = v;
		b = t;
		break;

	case 3:
		r = p;
		g = q;
		b = v;
		break;
	case 4:
		r = t;
		g = p;
		b = v;
		break;
	//     case 5:
	default:
		r = v;
		g = p;
		b = q;
		break;
	}
}

class path_indx_class
{
  public:
	int path_id;
	std::size_t indx;
	path_indx_class()
	{
	}
	inline bool operator<(const path_indx_class &rhs) { return path_id < rhs.path_id; }
	inline bool operator>(const path_indx_class &rhs) { return path_id > rhs.path_id; }
};

template <typename T>
void _mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	try
	{

		Vector<std::size_t, 3> target_shape;
		target_shape = std::size_t(0);
		const mxArray *params = prhs[nrhs - 1];
		T b = 5;
		bool consider_scale = false;
		int mode = 0;
		bool weight = true;
		bool normalize = true;
		T sigma = 0.25;
		T gamma = 0.5;
		Vector<T, 3> offset;
		Vector<T, 3> vscale;
		offset = T(0);
		vscale = T(1);

		bool lininterp = true;

		if (mxIsCell(params))
		{
			if (mhs::mex_hasParam(params, "shape") != -1)
			{
				target_shape = mhs::mex_getParam<T>(params, "shape", 3);
			}
			if (mhs::mex_hasParam(params, "offset") != -1)
			{
				offset = mhs::mex_getParam<T>(params, "offset", 3);
			}
			if (mhs::mex_hasParam(params, "vscale") != -1)
			{
				vscale = mhs::mex_getParam<T>(params, "vscale", 3);
			}
			if (mhs::mex_hasParam(params, "b") != -1)
			{
				b = mhs::mex_getParam<T>(params, "b", 1)[0];
			}

			if (mhs::mex_hasParam(params, "consider_scale") != -1)
			{
				consider_scale = mhs::mex_getParam<bool>(params, "consider_scale", 1)[0];
			}

			if (mhs::mex_hasParam(params, "mode") != -1)
			{
				mode = mhs::mex_getParam<int>(params, "mode", 1)[0];
			}
			if (mhs::mex_hasParam(params, "sigma") != -1)
			{
				sigma = mhs::mex_getParam<T>(params, "sigma", 1)[0];
			}
			if (mhs::mex_hasParam(params, "gamma") != -1)
			{
				gamma = mhs::mex_getParam<T>(params, "gamma", 1)[0];
			}

			if (mhs::mex_hasParam(params, "weight") != -1)
			{
				weight = mhs::mex_getParam<bool>(params, "weight", 1)[0];
			}

			if (mhs::mex_hasParam(params, "normalize") != -1)
			{
				normalize = mhs::mex_getParam<bool>(params, "normalize", 1)[0];
			}

			if (mhs::mex_hasParam(params, "lininterp") != -1)
			{
				lininterp = mhs::mex_getParam<bool>(params, "lininterp", 1)[0];
			}
		}
		//

		offset.print();

		mhs::dataArray<T> edges;
		edges = mhs::dataArray<T>(prhs[0]);

		T *result = NULL; //( T * ) mxGetData ( plhs[0]);
		T *result_rgb = NULL;
		T *result_gray = NULL;

		T *accu = NULL;
		int *accu_id = NULL;

		switch (mode)
		{
		case 0:
		{
			mwSize ndims[3];
			ndims[2] = target_shape[0];
			ndims[1] = target_shape[1];
			ndims[0] = target_shape[2];

			plhs[0] = mxCreateNumericArray(3, ndims, mhs::mex_getClassId<T>(), mxREAL);
			result = (T *)mxGetData(plhs[0]);
		}
		break;

		case 1:
		{
			mwSize ndims[3];
			ndims[2] = target_shape[0];
			ndims[1] = target_shape[1];
			ndims[0] = target_shape[2];

			plhs[0] = mxCreateNumericArray(3, ndims, mhs::mex_getClassId<T>(), mxREAL);
			result_gray = (T *)mxGetData(plhs[0]);

			accu = new T[ndims[2] * ndims[1] * ndims[0]];
			accu_id = new int[ndims[2] * ndims[1] * ndims[0]];
		}
		break;

		case 2:
		{
			mwSize ndims[4];
			ndims[3] = target_shape[0];
			ndims[2] = target_shape[1];
			ndims[1] = target_shape[2];
			ndims[0] = 3;

			accu = new T[ndims[3] * ndims[2] * ndims[1] * 4];
			accu_id = new int[ndims[3] * ndims[2] * ndims[1]];

			plhs[0] = mxCreateNumericArray(4, ndims, mhs::mex_getClassId<T>(), mxREAL);

			result_rgb = (T *)mxGetData(plhs[0]);
		}
		break;

		default:
			return;
		}

		sta_assert_error(edges.dim.size() == 2);
		sta_assert_error(edges.dim[0] % 2 == 0);

		std::size_t num_edges = edges.dim[0] / 2;

		std::size_t numv = target_shape[0] * target_shape[1] * target_shape[2];

		if (mode == 0)
		{
#pragma omp parallel for num_threads(omp_get_num_procs())
			for (std::size_t e = 0; e < numv; e++)
			{
				result[e] = b;
			}
		}
		if (mode > 0)
		{
			for (std::size_t e = 0; e < numv; e++)
			{
				accu_id[e] = -1000;
			}
		}

		const T *data = edges.data;

		int path_ids[2];
		path_ids[0] = std::numeric_limits<int>::max();
		path_ids[1] = -std::numeric_limits<int>::max();

		std::vector<path_indx_class> indeces;
		indeces.resize(num_edges);
		for (std::size_t e = 0; e < num_edges; e++)
		{
			const T *p_data = data + e * 10;
			int path_id = p_data[4];
			indeces[e].path_id = path_id;
			indeces[e].indx = e;
		}

		std::sort(indeces.begin(), indeces.end());

		if (mode > 0)
		{
			for (std::size_t e = 0; e < num_edges; e++)
			{
				const T *p_data = data + e * 10;
				int path_id = p_data[4];
				path_ids[0] = std::min(path_ids[0], path_id);
				path_ids[1] = std::max(path_ids[1], path_id);
			}
		}

#pragma omp parallel for num_threads(omp_get_num_procs())
		for (std::size_t e = 0; e < num_edges; e++)
		{

			const T *p_data = data + indeces[e].indx * 10;
			int path_id = indeces[e].path_id;

			Vector<T, 3> x1(p_data[2], p_data[1], p_data[0]);
			T s1 = p_data[3];
			Vector<T, 3> x2(p_data[7], p_data[6], p_data[5]);
			T s2 = p_data[8];

			x1 += offset;
			x2 += offset;

			int path_id2 = p_data[4];
			sta_assert_error(path_id2 == path_id);

			std::size_t bb[2][3];

			for (int i = 0; i < 3; i++)
			{
				int b2 = std::ceil(b * vscale[i]);
				if (!consider_scale)
				{
					bb[0][i] = std::max(std::min(x1[i] - b2, x2[i] - b2), T(0));
					bb[1][i] = std::max(std::min(std::max(x1[i] + b2, x2[i] + b2), T(target_shape[i])), T(0));
				}
				else
				{
					bb[0][i] = std::max(std::min(x1[i] - s1 - b2 * s1, x2[i] - s2 - b2 * s2), T(0));
					bb[1][i] = std::max(std::min(std::max(x1[i] + s1 + b2 * s1, x2[i] + s2 + b2 * s2), T(target_shape[i])), T(0));
				}
			}

			Vector<T, 3> x3;

			for (std::size_t z = bb[0][0]; z < bb[1][0]; z++)
			{
				x3[0] = z;
				for (std::size_t y = bb[0][1]; y < bb[1][1]; y++)
				{
					x3[1] = y;
					for (std::size_t x = bb[0][2]; x < bb[1][2]; x++)
					{
						x3[2] = x;
						if (true)
						{

							std::size_t indx = SUB2IND(x, y, z, target_shape);
							sta_assert_error(indx < numv);
							T score = 0;
							Vector<T, 3> x1_s = x1 / vscale;
							Vector<T, 3> x2_s = x2 / vscale;
							Vector<T, 3> x3_s = x3 / vscale;
							if (!consider_scale)
							{
								score = get_dist_from_line(x1_s, x2_s, x3_s);
							}
							else
							{
								score = get_dist_from_surface(x1, s1, x2, s2, x3);
							}
							switch (mode)
							{
							case 0:
							{
								sta_assert_error(result != NULL);
								result[indx] = std::min(result[indx], score);
							}
							break;
							case 1:
							{
								sta_assert_error(result_gray != NULL);
								float alpha;
								if (lininterp)
								{
									alpha = std::max(sigma - score, T(0));
								}
								else
								{
									alpha = std::exp(-score * score / (sigma * sigma));
								}

								if (weight)
								{
									if (accu_id[indx] == path_id)
									{
										if (accu[indx] < alpha)
										{
											result_gray[indx] += alpha - accu[indx];
											accu[indx] = alpha;
										}
									}
									else
									{
										result_gray[indx] += alpha;
										accu[indx] = alpha;
										accu_id[indx] = path_id;
									}
								}
								else
								{
									if (alpha > result_gray[indx])
									{
										result_gray[indx] = alpha;
									}
								}
							}
							break;

							case 2:
							{
								sta_assert_error(result_rgb != NULL);
								float alpha;
								if (lininterp)
								{
									alpha = std::max(sigma - score, T(0));
								}
								else
								{
									alpha = std::exp(-score * score / (sigma * sigma));
								}
								Vector<T, 3> color;
								color = x1 - x2;
								color.normalize();
								color[0] = std::abs(color[0]);
								color[1] = std::abs(color[1]);
								color[2] = std::abs(color[2]);

								{
									if (weight)
									{
										sta_assert_error(((indx * 4 + 3) < target_shape[0] * target_shape[1] * target_shape[2] * 4)) if (accu_id[indx] == path_id)
										{
											if (accu[indx * 4] < alpha)
											{
												for (int c = 0; c < 3; c++)
												{
													sta_assert_error((indx * 4 + c + 1) < target_shape[0] * target_shape[1] * target_shape[2] * 4)
														result_rgb[indx * 3 + c] += color[c] * alpha - accu[indx * 4 + c + 1];
													accu[indx * 4 + c + 1] = color[c] * alpha;
												}
												accu[indx * 4] = alpha;
											}
										}
										else
										{
											for (int c = 0; c < 3; c++)
											{
												result_rgb[indx * 3 + c] += color[c] * alpha;
												sta_assert_error((indx * 4 + c + 1) < target_shape[0] * target_shape[1] * target_shape[2] * 4)
													accu[indx * 4 + c + 1] = color[c] * alpha;
											}
											accu[indx * 4] = alpha;
											accu_id[indx] = path_id;
										}
									}
									else
									{
										T old_score = result_rgb[indx * 3 + 0] * result_rgb[indx * 3 + 0] +
													  result_rgb[indx * 3 + 1] * result_rgb[indx * 3 + 1] +
													  result_rgb[indx * 3 + 2] * result_rgb[indx * 3 + 2];
										if (alpha * alpha > old_score)
										{

											for (int c = 0; c < 3; c++)
												result_rgb[indx * 3 + c] = result_rgb[indx * 3 + c] * (1 - alpha) + color[c] * (alpha);
										}
									}
								}
							}
							break;
							}
						}
					}
				}
			}
		}

		if (weight && normalize)
		{
			T maxw = 0.000000000000001;
			T maxw_g = 0.000000000000001;

#pragma omp parallel for num_threads(omp_get_num_procs())
			for (std::size_t e = 0; e < numv; e++)
			{
				if (result_rgb != NULL)
				{
					T *c = result_rgb + e * 3;
					T n = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);

					T nsqrt = n;
					if (std::abs(gamma - 1) > 0.0000001)
					{
						nsqrt = std::pow(n, gamma);
					}
					T w = nsqrt / (n + 0.00000000000000000000001);
					c[0] *= w;
					c[1] *= w;
					c[2] *= w;
					maxw = std::max(maxw, nsqrt);
				}
				if (result_gray != NULL)
				{
					maxw_g = std::max(maxw_g, result_gray[e]);
				}
				//maxw=std::max(maxw,n);
			}

#pragma omp parallel for num_threads(omp_get_num_procs())
			for (std::size_t e = 0; e < numv; e++)
			{
				if (result_rgb != NULL)
				{
					T *c = result_rgb + e * 3;
					c[0] /= maxw;
					c[1] /= maxw;
					c[2] /= maxw;
				}
				if (result_gray != NULL)
				{
					result_gray[e] /= maxw_g;
				}
			}
		}

		if (accu != NULL)
		{
			delete[] accu;
			delete[] accu_id;
		}
	}
	catch (mhs::STAError &error)
	{
		printf("error cleaning up!!\n");
		mexErrMsgTxt(error.what());
	}
	catch (...)
	{
		printf("error cleaning up!!\n");
		mexErrMsgTxt("exeption ?");
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs < 2)
		mexErrMsgTxt("error: nrhs<2\n");

	const mxArray *img = prhs[0];

	if ((img == NULL) || (mxGetClassID(img) == mxDOUBLE_CLASS))
		_mexFunction<double>(nlhs, plhs, nrhs, prhs);
	else if (mxGetClassID(img) == mxSINGLE_CLASS)
		_mexFunction<float>(nlhs, plhs, nrhs, prhs);
	else
		mexErrMsgTxt("error: unsupported data type\n");
}
