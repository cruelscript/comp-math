#pragma once
#include <cmath>
#include "cmath.h"

#define SIGN0(x) (((x) > 0) ? 1 : -1)
#define SIGN(x) (((x) == 0) ? 0 : SIGN0(x))
#define SIGN2(a, b) (SIGN(b)*abs(a))

double FMin(double(*F)(double), double a, double b, double tol)
{
	double c = (double)((double)3.0 - (double)sqrt(5.0)) / (double)2,//с - это возведённая в квадрат величина, обратная к золотому сечению
		d, e = (double)0, eps = (double)sqrt(EPSILON),//eps приблизительно равно квадратному корню из относительной машинной точности
		xm, p, q, r, tol1, tol2, u, v = a + c * (b - a), w = v, x = v, fx = F(v), fu, fv = fx, fw = fx;

	while (1)
	{
		xm = (a + b) / (double)2;
		tol1 = eps * abs(x) + tol / (double)3;
		tol2 = tol1 * (double)2;

		//Проверить критерий окончания
		if (abs(x - xm) <= (tol2 - (b - a) / (double)2)) return x;

		if (abs(e) > tol1)//Если золотое сечение не требуется
		{
			//Построить параболу
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = (double)2 * (q - r);
			if (q > 0) p = -p;
			else q = -q;//q = ABS(q)
			r = e;
			e = d;

			//Приемлема ли парабола
			if (abs(p) >= abs(q * r / (double)2) || p <= q * (a - x) || p >= q * (b - x))
			{
				//Шаг золотого сечения
				if (x >= xm) e = a - x;
				else e = b - x;
				d = c * e;
			}
			else
			{
				//Шаг параболической интерполяции
				d = p / q;
				u = x + d;
				//F не следует вычислять слишком близко к 'a' или 'b'
				if ((u - a) < tol2 || (b - u) < tol2) d = SIGN2(tol1, xm - x);
			}
		}
		else
		{
			//Шаг золотого сечения
			if (x >= xm) e = a - x;
			else e = b - x;
			d = c * e;
		}

		//F не следует вычислять слишком близко к 'x'
		if (abs(d) >= tol1) u = x + d;
		else u = x + SIGN2(tol1, d);
		fu = F(u);

		//Присвоить новые значения параметрам 'a', 'b', 'v', 'w' и 'x'
		if (fu <= fx)
		{
			if (u >= x) a = x;
			else b = x;
			v = w;
			fv = fw;
			w = x;
			fw = fx;
			x = u;
			fx = fu;
			continue;
		}
		if (u < x) a = u;
		else b = u;
		if (fu <= fw || w == x)
		{
			v = w;
			fv = fw;
			w = u;
			fw = fu;
		}
		else if (fu <= fv || v == x || v == w)
		{
			v = u;
			fv = fu;
		}
	}

	return x;
}
