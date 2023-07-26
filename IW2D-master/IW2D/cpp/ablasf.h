/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#ifndef _ablasf_h
#define _ablasf_h

#include "ap.h"
#include "amp.h"
namespace ablasf
{
    template<unsigned int Precision>
    bool cmatrixrank1f(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::campf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int iv);
    template<unsigned int Precision>
    bool rmatrixrank1f(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::ampf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        int iv);
    template<unsigned int Precision>
    bool cmatrixmvf(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int opa,
        ap::template_1d_array< amp::campf<Precision> >& x,
        int ix,
        ap::template_1d_array< amp::campf<Precision> >& y,
        int iy);
    template<unsigned int Precision>
    bool rmatrixmvf(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int opa,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        int ix,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int iy);
    template<unsigned int Precision>
    bool cmatrixrighttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    bool cmatrixlefttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    bool rmatrixrighttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    bool rmatrixlefttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    bool cmatrixsyrkf(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc,
        bool isupper);
    template<unsigned int Precision>
    bool rmatrixsyrkf(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc,
        bool isupper);
    template<unsigned int Precision>
    bool rmatrixgemmf(int m,
        int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc);
    template<unsigned int Precision>
    bool cmatrixgemmf(int m,
        int n,
        int k,
        amp::campf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::campf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc);


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixrank1f(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::campf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int iv)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_cmatrixrank1f<Precision>(m, n, a, ia, ja, u, iu, v, iv);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixrank1f(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::ampf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        int iv)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_rmatrixrank1f<Precision>(m, n, a, ia, ja, u, iu, v, iv);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixmvf(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int opa,
        ap::template_1d_array< amp::campf<Precision> >& x,
        int ix,
        ap::template_1d_array< amp::campf<Precision> >& y,
        int iy)
    {
        bool result;


        result = false;
        return result;
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixmvf(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int opa,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        int ix,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int iy)
    {
        bool result;


        result = false;
        return result;
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixrighttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_cmatrixrighttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixlefttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_cmatrixlefttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixrighttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_rmatrixrighttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixlefttrsmf(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_rmatrixlefttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixsyrkf(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc,
        bool isupper)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_cmatrixsyrkf<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixsyrkf(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc,
        bool isupper)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_rmatrixsyrkf<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixgemmf(int m,
        int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_rmatrixgemmf<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
    #endif
    }


    /*************************************************************************
    Fast kernel

      -- ALGLIB routine --
         19.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixgemmf(int m,
        int n,
        int k,
        amp::campf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::campf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_ABLAS
        bool result;


        result = false;
        return result;
    #else
        return amp::_i_cmatrixgemmf<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
    #endif
    }
} // namespace

#endif
