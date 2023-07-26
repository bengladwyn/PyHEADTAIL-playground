/*************************************************************************
Copyright (c) 2009-2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _ablas_h
#define _ablas_h

#include "ap.h"
#include "amp.h"
#include "ablasf.h"
namespace ablas
{
    template<unsigned int Precision>
    void ablassplitlength(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int& n1,
        int& n2);
    template<unsigned int Precision>
    void ablascomplexsplitlength(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        int& n1,
        int& n2);
    template<unsigned int Precision>
    int ablasblocksize(const ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    int ablascomplexblocksize(const ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    int ablasmicroblocksize();
    template<unsigned int Precision>
    void cmatrixtranspose(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb);
    template<unsigned int Precision>
    void rmatrixtranspose(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb);
    template<unsigned int Precision>
    void cmatrixcopy(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb);
    template<unsigned int Precision>
    void rmatrixcopy(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb);
    template<unsigned int Precision>
    void cmatrixrank1(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::campf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int iv);
    template<unsigned int Precision>
    void rmatrixrank1(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::ampf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        int iv);
    template<unsigned int Precision>
    void cmatrixmv(int m,
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
    void rmatrixmv(int m,
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
    void cmatrixrighttrsm(int m,
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
    void cmatrixlefttrsm(int m,
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
    void rmatrixrighttrsm(int m,
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
    void rmatrixlefttrsm(int m,
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
    void cmatrixsyrk(int n,
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
    void rmatrixsyrk(int n,
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
    void cmatrixgemm(int m,
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
    template<unsigned int Precision>
    void rmatrixgemm(int m,
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
    void ablasinternalsplitlength(int n,
        int nb,
        int& n1,
        int& n2);
    template<unsigned int Precision>
    void cmatrixrighttrsm2(int m,
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
    void cmatrixlefttrsm2(int m,
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
    void rmatrixrighttrsm2(int m,
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
    void rmatrixlefttrsm2(int m,
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
    void cmatrixsyrk2(int n,
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
    void rmatrixsyrk2(int n,
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
    void cmatrixgemmk(int m,
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
    template<unsigned int Precision>
    void rmatrixgemmk(int m,
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


    /*************************************************************************
    Splits matrix length in two parts, left part should match ABLAS block size

    INPUT PARAMETERS
        A   -   real matrix, is passed to ensure that we didn't split
                complex matrix using real splitting subroutine.
                matrix itself is not changed.
        N   -   length, N>0

    OUTPUT PARAMETERS
        N1  -   length
        N2  -   length

    N1+N2=N, N1>=N2, N2 may be zero

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ablassplitlength(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int& n1,
        int& n2)
    {
        if( n>ablasblocksize<Precision>(a) )
        {
            ablasinternalsplitlength<Precision>(n, ablasblocksize<Precision>(a), n1, n2);
        }
        else
        {
            ablasinternalsplitlength<Precision>(n, ablasmicroblocksize<Precision>(), n1, n2);
        }
    }


    /*************************************************************************
    Complex ABLASSplitLength

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ablascomplexsplitlength(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        int& n1,
        int& n2)
    {
        if( n>ablascomplexblocksize<Precision>(a) )
        {
            ablasinternalsplitlength<Precision>(n, ablascomplexblocksize<Precision>(a), n1, n2);
        }
        else
        {
            ablasinternalsplitlength<Precision>(n, ablasmicroblocksize<Precision>(), n1, n2);
        }
    }


    /*************************************************************************
    Returns block size - subdivision size where  cache-oblivious  soubroutines
    switch to the optimized kernel.

    INPUT PARAMETERS
        A   -   real matrix, is passed to ensure that we didn't split
                complex matrix using real splitting subroutine.
                matrix itself is not changed.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int ablasblocksize(const ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        int result;


        result = 32;
        return result;
    }


    /*************************************************************************
    Block size for complex subroutines.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int ablascomplexblocksize(const ap::template_2d_array< amp::campf<Precision> >& a)
    {
        int result;


        result = 24;
        return result;
    }


    /*************************************************************************
    Microblock size

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int ablasmicroblocksize()
    {
        int result;


        result = 8;
        return result;
    }


    /*************************************************************************
    Cache-oblivous complex "copy-and-transpose"

    Input parameters:
        M   -   number of rows
        N   -   number of columns
        A   -   source matrix, MxN submatrix is copied and transposed
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        A   -   destination matrix
        IB  -   submatrix offset (row index)
        JB  -   submatrix offset (column index)
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixtranspose(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb)
    {
        int i;
        int s1;
        int s2;
        int i_;
        int i1_;


        if( m<=2*ablascomplexblocksize<Precision>(a) && n<=2*ablascomplexblocksize<Precision>(a) )
        {
            
            //
            // base case
            //
            for(i=0; i<=m-1; i++)
            {
                i1_ = (ja) - (ib);
                for(i_=ib; i_<=ib+n-1;i_++)
                {
                    b(i_,jb+i) = a(ia+i,i_+i1_);
                }
            }
        }
        else
        {
            
            //
            // Cache-oblivious recursion
            //
            if( m>n )
            {
                ablascomplexsplitlength<Precision>(a, m, s1, s2);
                cmatrixtranspose<Precision>(s1, n, a, ia, ja, b, ib, jb);
                cmatrixtranspose<Precision>(s2, n, a, ia+s1, ja, b, ib, jb+s1);
            }
            else
            {
                ablascomplexsplitlength<Precision>(a, n, s1, s2);
                cmatrixtranspose<Precision>(m, s1, a, ia, ja, b, ib, jb);
                cmatrixtranspose<Precision>(m, s2, a, ia, ja+s1, b, ib+s1, jb);
            }
        }
    }


    /*************************************************************************
    Cache-oblivous real "copy-and-transpose"

    Input parameters:
        M   -   number of rows
        N   -   number of columns
        A   -   source matrix, MxN submatrix is copied and transposed
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        A   -   destination matrix
        IB  -   submatrix offset (row index)
        JB  -   submatrix offset (column index)
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixtranspose(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb)
    {
        int i;
        int s1;
        int s2;


        if( m<=2*ablasblocksize<Precision>(a) && n<=2*ablasblocksize<Precision>(a) )
        {
            
            //
            // base case
            //
            for(i=0; i<=m-1; i++)
            {
                amp::vmove(b.getcolumn(jb+i, ib, ib+n-1), a.getrow(ia+i, ja, ja+n-1));
            }
        }
        else
        {
            
            //
            // Cache-oblivious recursion
            //
            if( m>n )
            {
                ablassplitlength<Precision>(a, m, s1, s2);
                rmatrixtranspose<Precision>(s1, n, a, ia, ja, b, ib, jb);
                rmatrixtranspose<Precision>(s2, n, a, ia+s1, ja, b, ib, jb+s1);
            }
            else
            {
                ablassplitlength<Precision>(a, n, s1, s2);
                rmatrixtranspose<Precision>(m, s1, a, ia, ja, b, ib, jb);
                rmatrixtranspose<Precision>(m, s2, a, ia, ja+s1, b, ib+s1, jb);
            }
        }
    }


    /*************************************************************************
    Copy

    Input parameters:
        M   -   number of rows
        N   -   number of columns
        A   -   source matrix, MxN submatrix is copied and transposed
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        B   -   destination matrix
        IB  -   submatrix offset (row index)
        JB  -   submatrix offset (column index)
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixcopy(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb)
    {
        int i;
        int i_;
        int i1_;


        for(i=0; i<=m-1; i++)
        {
            i1_ = (ja) - (jb);
            for(i_=jb; i_<=jb+n-1;i_++)
            {
                b(ib+i,i_) = a(ia+i,i_+i1_);
            }
        }
    }


    /*************************************************************************
    Copy

    Input parameters:
        M   -   number of rows
        N   -   number of columns
        A   -   source matrix, MxN submatrix is copied and transposed
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        B   -   destination matrix
        IB  -   submatrix offset (row index)
        JB  -   submatrix offset (column index)
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixcopy(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb)
    {
        int i;


        for(i=0; i<=m-1; i++)
        {
            amp::vmove(b.getrow(ib+i, jb, jb+n-1), a.getrow(ia+i, ja, ja+n-1));
        }
    }


    /*************************************************************************
    Rank-1 correction: A := A + u*v'

    INPUT PARAMETERS:
        M   -   number of rows
        N   -   number of columns
        A   -   target matrix, MxN submatrix is updated
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        U   -   vector #1
        IU  -   subvector offset
        V   -   vector #2
        IV  -   subvector offset
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrank1(int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::campf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int iv)
    {
        int i;
        amp::campf<Precision> s;
        int i_;
        int i1_;


        if( m==0 || n==0 )
        {
            return;
        }
        if( ablasf::cmatrixrank1f<Precision>(m, n, a, ia, ja, u, iu, v, iv) )
        {
            return;
        }
        for(i=0; i<=m-1; i++)
        {
            s = u(iu+i);
            i1_ = (iv) - (ja);
            for(i_=ja; i_<=ja+n-1;i_++)
            {
                a(ia+i,i_) = a(ia+i,i_) + s*v(i_+i1_);
            }
        }
    }


    /*************************************************************************
    Rank-1 correction: A := A + u*v'

    INPUT PARAMETERS:
        M   -   number of rows
        N   -   number of columns
        A   -   target matrix, MxN submatrix is updated
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        U   -   vector #1
        IU  -   subvector offset
        V   -   vector #2
        IV  -   subvector offset
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrank1(int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        ap::template_1d_array< amp::ampf<Precision> >& u,
        int iu,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        int iv)
    {
        int i;
        amp::ampf<Precision> s;


        if( m==0 || n==0 )
        {
            return;
        }
        if( ablasf::rmatrixrank1f<Precision>(m, n, a, ia, ja, u, iu, v, iv) )
        {
            return;
        }
        for(i=0; i<=m-1; i++)
        {
            s = u(iu+i);
            amp::vadd(a.getrow(ia+i, ja, ja+n-1), v.getvector(iv, iv+n-1), s);
        }
    }


    /*************************************************************************
    Matrix-vector product: y := op(A)*x

    INPUT PARAMETERS:
        M   -   number of rows of op(A)
                M>=0
        N   -   number of columns of op(A)
                N>=0
        A   -   target matrix
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        OpA -   operation type:
                * OpA=0     =>  op(A) = A
                * OpA=1     =>  op(A) = A^T
                * OpA=2     =>  op(A) = A^H
        X   -   input vector
        IX  -   subvector offset
        IY  -   subvector offset

    OUTPUT PARAMETERS:
        Y   -   vector which stores result

    if M=0, then subroutine does nothing.
    if N=0, Y is filled by zeros.


      -- ALGLIB routine --

         28.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixmv(int m,
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
        int i;
        amp::campf<Precision> v;
        int i_;
        int i1_;


        if( m==0 )
        {
            return;
        }
        if( n==0 )
        {
            for(i=0; i<=m-1; i++)
            {
                y(iy+i) = 0;
            }
            return;
        }
        if( ablasf::cmatrixmvf<Precision>(m, n, a, ia, ja, opa, x, ix, y, iy) )
        {
            return;
        }
        if( opa==0 )
        {
            
            //
            // y = A*x
            //
            for(i=0; i<=m-1; i++)
            {
                i1_ = (ix)-(ja);
                v = 0.0;
                for(i_=ja; i_<=ja+n-1;i_++)
                {
                    v += a(ia+i,i_)*x(i_+i1_);
                }
                y(iy+i) = v;
            }
            return;
        }
        if( opa==1 )
        {
            
            //
            // y = A^T*x
            //
            for(i=0; i<=m-1; i++)
            {
                y(iy+i) = 0;
            }
            for(i=0; i<=n-1; i++)
            {
                v = x(ix+i);
                i1_ = (ja) - (iy);
                for(i_=iy; i_<=iy+m-1;i_++)
                {
                    y(i_) = y(i_) + v*a(ia+i,i_+i1_);
                }
            }
            return;
        }
        if( opa==2 )
        {
            
            //
            // y = A^H*x
            //
            for(i=0; i<=m-1; i++)
            {
                y(iy+i) = 0;
            }
            for(i=0; i<=n-1; i++)
            {
                v = x(ix+i);
                i1_ = (ja) - (iy);
                for(i_=iy; i_<=iy+m-1;i_++)
                {
                    y(i_) = y(i_) + v*amp::conj(a(ia+i,i_+i1_));
                }
            }
            return;
        }
    }


    /*************************************************************************
    Matrix-vector product: y := op(A)*x

    INPUT PARAMETERS:
        M   -   number of rows of op(A)
        N   -   number of columns of op(A)
        A   -   target matrix
        IA  -   submatrix offset (row index)
        JA  -   submatrix offset (column index)
        OpA -   operation type:
                * OpA=0     =>  op(A) = A
                * OpA=1     =>  op(A) = A^T
        X   -   input vector
        IX  -   subvector offset
        IY  -   subvector offset

    OUTPUT PARAMETERS:
        Y   -   vector which stores result

    if M=0, then subroutine does nothing.
    if N=0, Y is filled by zeros.


      -- ALGLIB routine --

         28.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixmv(int m,
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
        int i;
        amp::ampf<Precision> v;


        if( m==0 )
        {
            return;
        }
        if( n==0 )
        {
            for(i=0; i<=m-1; i++)
            {
                y(iy+i) = 0;
            }
            return;
        }
        if( ablasf::rmatrixmvf<Precision>(m, n, a, ia, ja, opa, x, ix, y, iy) )
        {
            return;
        }
        if( opa==0 )
        {
            
            //
            // y = A*x
            //
            for(i=0; i<=m-1; i++)
            {
                v = amp::vdotproduct(a.getrow(ia+i, ja, ja+n-1), x.getvector(ix, ix+n-1));
                y(iy+i) = v;
            }
            return;
        }
        if( opa==1 )
        {
            
            //
            // y = A^T*x
            //
            for(i=0; i<=m-1; i++)
            {
                y(iy+i) = 0;
            }
            for(i=0; i<=n-1; i++)
            {
                v = x(ix+i);
                amp::vadd(y.getvector(iy, iy+m-1), a.getrow(ia+i, ja, ja+m-1), v);
            }
            return;
        }
    }


    /*************************************************************************
    This subroutine calculates X*op(A^-1) where:
    * X is MxN general matrix
    * A is NxN upper/lower triangular/unitriangular matrix
    * "op" may be identity transformation, transposition, conjugate transposition

    Multiplication result replaces X.
    Cache-oblivious algorithm is used.

    INPUT PARAMETERS
        N   -   matrix size, N>=0
        M   -   matrix size, N>=0
        A       -   matrix, actial matrix is stored in A[I1:I1+N-1,J1:J1+N-1]
        I1      -   submatrix offset
        J1      -   submatrix offset
        IsUpper -   whether matrix is upper triangular
        IsUnit  -   whether matrix is unitriangular
        OpType  -   transformation type:
                    * 0 - no transformation
                    * 1 - transposition
                    * 2 - conjugate transposition
        C   -   matrix, actial matrix is stored in C[I2:I2+M-1,J2:J2+N-1]
        I2  -   submatrix offset
        J2  -   submatrix offset

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrighttrsm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablascomplexblocksize<Precision>(a);
        if( m<=bs && n<=bs )
        {
            cmatrixrighttrsm2<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( m>=n )
        {
            
            //
            // Split X: X*A = (X1 X2)^T*A
            //
            ablascomplexsplitlength<Precision>(a, m, s1, s2);
            cmatrixrighttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixrighttrsm<Precision>(s2, n, a, i1, j1, isupper, isunit, optype, x, i2+s1, j2);
        }
        else
        {
            
            //
            // Split A:
            //               (A1  A12)
            // X*op(A) = X*op(       )
            //               (     A2)
            //
            // Different variants depending on
            // IsUpper/OpType combinations
            //
            ablascomplexsplitlength<Precision>(a, n, s1, s2);
            if( isupper && optype==0 )
            {
                
                //
                //                  (A1  A12)-1
                // X*A^-1 = (X1 X2)*(       )
                //                  (     A2)
                //
                cmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                cmatrixgemm<Precision>(m, s2, s1, -amp::ampf<Precision>("1.0"), x, i2, j2, 0, a, i1, j1+s1, 0, amp::ampf<Precision>("1.0"), x, i2, j2+s1);
                cmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                return;
            }
            if( isupper && optype!=0 )
            {
                
                //
                //                  (A1'     )-1
                // X*A^-1 = (X1 X2)*(        )
                //                  (A12' A2')
                //
                cmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                cmatrixgemm<Precision>(m, s1, s2, -amp::ampf<Precision>("1.0"), x, i2, j2+s1, 0, a, i1, j1+s1, optype, amp::ampf<Precision>("1.0"), x, i2, j2);
                cmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( !isupper && optype==0 )
            {
                
                //
                //                  (A1     )-1
                // X*A^-1 = (X1 X2)*(       )
                //                  (A21  A2)
                //
                cmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                cmatrixgemm<Precision>(m, s1, s2, -amp::ampf<Precision>("1.0"), x, i2, j2+s1, 0, a, i1+s1, j1, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                cmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( !isupper && optype!=0 )
            {
                
                //
                //                  (A1' A21')-1
                // X*A^-1 = (X1 X2)*(        )
                //                  (     A2')
                //
                cmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                cmatrixgemm<Precision>(m, s2, s1, -amp::ampf<Precision>("1.0"), x, i2, j2, 0, a, i1+s1, j1, optype, amp::ampf<Precision>("1.0"), x, i2, j2+s1);
                cmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                return;
            }
        }
    }


    /*************************************************************************
    This subroutine calculates op(A^-1)*X where:
    * X is MxN general matrix
    * A is MxM upper/lower triangular/unitriangular matrix
    * "op" may be identity transformation, transposition, conjugate transposition

    Multiplication result replaces X.
    Cache-oblivious algorithm is used.

    INPUT PARAMETERS
        N   -   matrix size, N>=0
        M   -   matrix size, N>=0
        A       -   matrix, actial matrix is stored in A[I1:I1+M-1,J1:J1+M-1]
        I1      -   submatrix offset
        J1      -   submatrix offset
        IsUpper -   whether matrix is upper triangular
        IsUnit  -   whether matrix is unitriangular
        OpType  -   transformation type:
                    * 0 - no transformation
                    * 1 - transposition
                    * 2 - conjugate transposition
        C   -   matrix, actial matrix is stored in C[I2:I2+M-1,J2:J2+N-1]
        I2  -   submatrix offset
        J2  -   submatrix offset

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlefttrsm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablascomplexblocksize<Precision>(a);
        if( m<=bs && n<=bs )
        {
            cmatrixlefttrsm2<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( n>=m )
        {
            
            //
            // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
            //
            ablascomplexsplitlength<Precision>(x, n, s1, s2);
            cmatrixlefttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixlefttrsm<Precision>(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2+s1);
        }
        else
        {
            
            //
            // Split A
            //
            ablascomplexsplitlength<Precision>(a, m, s1, s2);
            if( isupper && optype==0 )
            {
                
                //
                //           (A1  A12)-1  ( X1 )
                // A^-1*X* = (       )   *(    )
                //           (     A2)    ( X2 )
                //
                cmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                cmatrixgemm<Precision>(s1, n, s2, -amp::ampf<Precision>("1.0"), a, i1, j1+s1, 0, x, i2+s1, j2, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                cmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( isupper && optype!=0 )
            {
                
                //
                //          (A1'     )-1 ( X1 )
                // A^-1*X = (        )  *(    )
                //          (A12' A2')   ( X2 )
                //
                cmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                cmatrixgemm<Precision>(s2, n, s1, -amp::ampf<Precision>("1.0"), a, i1, j1+s1, optype, x, i2, j2, 0, amp::ampf<Precision>("1.0"), x, i2+s1, j2);
                cmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                return;
            }
            if( !isupper && optype==0 )
            {
                
                //
                //          (A1     )-1 ( X1 )
                // A^-1*X = (       )  *(    )
                //          (A21  A2)   ( X2 )
                //
                cmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                cmatrixgemm<Precision>(s2, n, s1, -amp::ampf<Precision>("1.0"), a, i1+s1, j1, 0, x, i2, j2, 0, amp::ampf<Precision>("1.0"), x, i2+s1, j2);
                cmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                return;
            }
            if( !isupper && optype!=0 )
            {
                
                //
                //          (A1' A21')-1 ( X1 )
                // A^-1*X = (        )  *(    )
                //          (     A2')   ( X2 )
                //
                cmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                cmatrixgemm<Precision>(s1, n, s2, -amp::ampf<Precision>("1.0"), a, i1+s1, j1, optype, x, i2+s1, j2, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                cmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
        }
    }


    /*************************************************************************
    Same as CMatrixRightTRSM, but for real matrices

    OpType may be only 0 or 1.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrighttrsm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablasblocksize<Precision>(a);
        if( m<=bs && n<=bs )
        {
            rmatrixrighttrsm2<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( m>=n )
        {
            
            //
            // Split X: X*A = (X1 X2)^T*A
            //
            ablassplitlength<Precision>(a, m, s1, s2);
            rmatrixrighttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixrighttrsm<Precision>(s2, n, a, i1, j1, isupper, isunit, optype, x, i2+s1, j2);
        }
        else
        {
            
            //
            // Split A:
            //               (A1  A12)
            // X*op(A) = X*op(       )
            //               (     A2)
            //
            // Different variants depending on
            // IsUpper/OpType combinations
            //
            ablassplitlength<Precision>(a, n, s1, s2);
            if( isupper && optype==0 )
            {
                
                //
                //                  (A1  A12)-1
                // X*A^-1 = (X1 X2)*(       )
                //                  (     A2)
                //
                rmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                rmatrixgemm<Precision>(m, s2, s1, -amp::ampf<Precision>("1.0"), x, i2, j2, 0, a, i1, j1+s1, 0, amp::ampf<Precision>("1.0"), x, i2, j2+s1);
                rmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                return;
            }
            if( isupper && optype!=0 )
            {
                
                //
                //                  (A1'     )-1
                // X*A^-1 = (X1 X2)*(        )
                //                  (A12' A2')
                //
                rmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                rmatrixgemm<Precision>(m, s1, s2, -amp::ampf<Precision>("1.0"), x, i2, j2+s1, 0, a, i1, j1+s1, optype, amp::ampf<Precision>("1.0"), x, i2, j2);
                rmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( !isupper && optype==0 )
            {
                
                //
                //                  (A1     )-1
                // X*A^-1 = (X1 X2)*(       )
                //                  (A21  A2)
                //
                rmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                rmatrixgemm<Precision>(m, s1, s2, -amp::ampf<Precision>("1.0"), x, i2, j2+s1, 0, a, i1+s1, j1, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                rmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( !isupper && optype!=0 )
            {
                
                //
                //                  (A1' A21')-1
                // X*A^-1 = (X1 X2)*(        )
                //                  (     A2')
                //
                rmatrixrighttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                rmatrixgemm<Precision>(m, s2, s1, -amp::ampf<Precision>("1.0"), x, i2, j2, 0, a, i1+s1, j1, optype, amp::ampf<Precision>("1.0"), x, i2, j2+s1);
                rmatrixrighttrsm<Precision>(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
                return;
            }
        }
    }


    /*************************************************************************
    Same as CMatrixLeftTRSM, but for real matrices

    OpType may be only 0 or 1.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlefttrsm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablasblocksize<Precision>(a);
        if( m<=bs && n<=bs )
        {
            rmatrixlefttrsm2<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( n>=m )
        {
            
            //
            // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
            //
            ablassplitlength<Precision>(x, n, s1, s2);
            rmatrixlefttrsm<Precision>(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixlefttrsm<Precision>(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2+s1);
        }
        else
        {
            
            //
            // Split A
            //
            ablassplitlength<Precision>(a, m, s1, s2);
            if( isupper && optype==0 )
            {
                
                //
                //           (A1  A12)-1  ( X1 )
                // A^-1*X* = (       )   *(    )
                //           (     A2)    ( X2 )
                //
                rmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                rmatrixgemm<Precision>(s1, n, s2, -amp::ampf<Precision>("1.0"), a, i1, j1+s1, 0, x, i2+s1, j2, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                rmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
            if( isupper && optype!=0 )
            {
                
                //
                //          (A1'     )-1 ( X1 )
                // A^-1*X = (        )  *(    )
                //          (A12' A2')   ( X2 )
                //
                rmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                rmatrixgemm<Precision>(s2, n, s1, -amp::ampf<Precision>("1.0"), a, i1, j1+s1, optype, x, i2, j2, 0, amp::ampf<Precision>("1.0"), x, i2+s1, j2);
                rmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                return;
            }
            if( !isupper && optype==0 )
            {
                
                //
                //          (A1     )-1 ( X1 )
                // A^-1*X = (       )  *(    )
                //          (A21  A2)   ( X2 )
                //
                rmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                rmatrixgemm<Precision>(s2, n, s1, -amp::ampf<Precision>("1.0"), a, i1+s1, j1, 0, x, i2, j2, 0, amp::ampf<Precision>("1.0"), x, i2+s1, j2);
                rmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                return;
            }
            if( !isupper && optype!=0 )
            {
                
                //
                //          (A1' A21')-1 ( X1 )
                // A^-1*X = (        )  *(    )
                //          (     A2')   ( X2 )
                //
                rmatrixlefttrsm<Precision>(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
                rmatrixgemm<Precision>(s1, n, s2, -amp::ampf<Precision>("1.0"), a, i1+s1, j1, optype, x, i2+s1, j2, 0, amp::ampf<Precision>("1.0"), x, i2, j2);
                rmatrixlefttrsm<Precision>(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
                return;
            }
        }
    }


    /*************************************************************************
    This subroutine calculates  C=alpha*A*A^H+beta*C  or  C=alpha*A^H*A+beta*C
    where:
    * C is NxN Hermitian matrix given by its upper/lower triangle
    * A is NxK matrix when A*A^H is calculated, KxN matrix otherwise

    Additional info:
    * cache-oblivious algorithm is used.
    * multiplication result replaces C. If Beta=0, C elements are not used in
      calculations (not multiplied by zero - just not referenced)
    * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
    * if both Beta and Alpha are zero, C is filled by zeros.

    INPUT PARAMETERS
        N       -   matrix size, N>=0
        K       -   matrix size, K>=0
        Alpha   -   coefficient
        A       -   matrix
        IA      -   submatrix offset
        JA      -   submatrix offset
        OpTypeA -   multiplication type:
                    * 0 - A*A^H is calculated
                    * 2 - A^H*A is calculated
        Beta    -   coefficient
        C       -   matrix
        IC      -   submatrix offset
        JC      -   submatrix offset
        IsUpper -   whether C is upper triangular or lower triangular

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixsyrk(int n,
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
        int s1;
        int s2;
        int bs;


        bs = ablascomplexblocksize<Precision>(a);
        if( n<=bs && k<=bs )
        {
            cmatrixsyrk2<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            return;
        }
        if( k>=n )
        {
            
            //
            // Split K
            //
            ablascomplexsplitlength<Precision>(a, k, s1, s2);
            if( optypea==0 )
            {
                cmatrixsyrk<Precision>(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixsyrk<Precision>(n, s2, alpha, a, ia, ja+s1, optypea, amp::ampf<Precision>("1.0"), c, ic, jc, isupper);
            }
            else
            {
                cmatrixsyrk<Precision>(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixsyrk<Precision>(n, s2, alpha, a, ia+s1, ja, optypea, amp::ampf<Precision>("1.0"), c, ic, jc, isupper);
            }
        }
        else
        {
            
            //
            // Split N
            //
            ablascomplexsplitlength<Precision>(a, n, s1, s2);
            if( optypea==0 && isupper )
            {
                cmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixgemm<Precision>(s1, s2, k, alpha, a, ia, ja, 0, a, ia+s1, ja, 2, beta, c, ic, jc+s1);
                cmatrixsyrk<Precision>(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea==0 && !isupper )
            {
                cmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixgemm<Precision>(s2, s1, k, alpha, a, ia+s1, ja, 0, a, ia, ja, 2, beta, c, ic+s1, jc);
                cmatrixsyrk<Precision>(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea!=0 && isupper )
            {
                cmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixgemm<Precision>(s1, s2, k, alpha, a, ia, ja, 2, a, ia, ja+s1, 0, beta, c, ic, jc+s1);
                cmatrixsyrk<Precision>(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea!=0 && !isupper )
            {
                cmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                cmatrixgemm<Precision>(s2, s1, k, alpha, a, ia, ja+s1, 2, a, ia, ja, 0, beta, c, ic+s1, jc);
                cmatrixsyrk<Precision>(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
        }
    }


    /*************************************************************************
    Same as CMatrixSYRK, but for real matrices

    OpType may be only 0 or 1.

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsyrk(int n,
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
        int s1;
        int s2;
        int bs;


        bs = ablasblocksize<Precision>(a);
        if( n<=bs && k<=bs )
        {
            rmatrixsyrk2<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            return;
        }
        if( k>=n )
        {
            
            //
            // Split K
            //
            ablassplitlength<Precision>(a, k, s1, s2);
            if( optypea==0 )
            {
                rmatrixsyrk<Precision>(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixsyrk<Precision>(n, s2, alpha, a, ia, ja+s1, optypea, amp::ampf<Precision>("1.0"), c, ic, jc, isupper);
            }
            else
            {
                rmatrixsyrk<Precision>(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixsyrk<Precision>(n, s2, alpha, a, ia+s1, ja, optypea, amp::ampf<Precision>("1.0"), c, ic, jc, isupper);
            }
        }
        else
        {
            
            //
            // Split N
            //
            ablassplitlength<Precision>(a, n, s1, s2);
            if( optypea==0 && isupper )
            {
                rmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixgemm<Precision>(s1, s2, k, alpha, a, ia, ja, 0, a, ia+s1, ja, 1, beta, c, ic, jc+s1);
                rmatrixsyrk<Precision>(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea==0 && !isupper )
            {
                rmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixgemm<Precision>(s2, s1, k, alpha, a, ia+s1, ja, 0, a, ia, ja, 1, beta, c, ic+s1, jc);
                rmatrixsyrk<Precision>(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea!=0 && isupper )
            {
                rmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixgemm<Precision>(s1, s2, k, alpha, a, ia, ja, 1, a, ia, ja+s1, 0, beta, c, ic, jc+s1);
                rmatrixsyrk<Precision>(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
            if( optypea!=0 && !isupper )
            {
                rmatrixsyrk<Precision>(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
                rmatrixgemm<Precision>(s2, s1, k, alpha, a, ia, ja+s1, 1, a, ia, ja, 0, beta, c, ic+s1, jc);
                rmatrixsyrk<Precision>(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
                return;
            }
        }
    }


    /*************************************************************************
    This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
    * C is MxN general matrix
    * op1(A) is MxK matrix
    * op2(B) is KxN matrix
    * "op" may be identity transformation, transposition, conjugate transposition

    Additional info:
    * cache-oblivious algorithm is used.
    * multiplication result replaces C. If Beta=0, C elements are not used in
      calculations (not multiplied by zero - just not referenced)
    * if Alpha=0, A is not used (not multiplied by zero - just not referenced)
    * if both Beta and Alpha are zero, C is filled by zeros.

    INPUT PARAMETERS
        N       -   matrix size, N>0
        M       -   matrix size, N>0
        K       -   matrix size, K>0
        Alpha   -   coefficient
        A       -   matrix
        IA      -   submatrix offset
        JA      -   submatrix offset
        OpTypeA -   transformation type:
                    * 0 - no transformation
                    * 1 - transposition
                    * 2 - conjugate transposition
        B       -   matrix
        IB      -   submatrix offset
        JB      -   submatrix offset
        OpTypeB -   transformation type:
                    * 0 - no transformation
                    * 1 - transposition
                    * 2 - conjugate transposition
        Beta    -   coefficient
        C       -   matrix
        IC      -   submatrix offset
        JC      -   submatrix offset

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixgemm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablascomplexblocksize<Precision>(a);
        if( m<=bs && n<=bs && k<=bs )
        {
            cmatrixgemmk<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            return;
        }
        if( m>=n && m>=k )
        {
            
            //
            // A*B = (A1 A2)^T*B
            //
            ablascomplexsplitlength<Precision>(a, m, s1, s2);
            cmatrixgemm<Precision>(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            if( optypea==0 )
            {
                cmatrixgemm<Precision>(s2, n, k, alpha, a, ia+s1, ja, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
            }
            else
            {
                cmatrixgemm<Precision>(s2, n, k, alpha, a, ia, ja+s1, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
            }
            return;
        }
        if( n>=m && n>=k )
        {
            
            //
            // A*B = A*(B1 B2)
            //
            ablascomplexsplitlength<Precision>(a, n, s1, s2);
            if( optypeb==0 )
            {
                cmatrixgemm<Precision>(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, s2, k, alpha, a, ia, ja, optypea, b, ib, jb+s1, optypeb, beta, c, ic, jc+s1);
            }
            else
            {
                cmatrixgemm<Precision>(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, s2, k, alpha, a, ia, ja, optypea, b, ib+s1, jb, optypeb, beta, c, ic, jc+s1);
            }
            return;
        }
        if( k>=m && k>=n )
        {
            
            //
            // A*B = (A1 A2)*(B1 B2)^T
            //
            ablascomplexsplitlength<Precision>(a, k, s1, s2);
            if( optypea==0 && optypeb==0 )
            {
                cmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib+s1, jb, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea==0 && optypeb!=0 )
            {
                cmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib, jb+s1, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea!=0 && optypeb==0 )
            {
                cmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib+s1, jb, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea!=0 && optypeb!=0 )
            {
                cmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                cmatrixgemm<Precision>(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib, jb+s1, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            return;
        }
    }


    /*************************************************************************
    Same as CMatrixGEMM, but for real numbers.
    OpType may be only 0 or 1.

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixgemm(int m,
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
        int s1;
        int s2;
        int bs;


        bs = ablasblocksize<Precision>(a);
        if( m<=bs && n<=bs && k<=bs )
        {
            rmatrixgemmk<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            return;
        }
        if( m>=n && m>=k )
        {
            
            //
            // A*B = (A1 A2)^T*B
            //
            ablassplitlength<Precision>(a, m, s1, s2);
            if( optypea==0 )
            {
                rmatrixgemm<Precision>(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(s2, n, k, alpha, a, ia+s1, ja, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
            }
            else
            {
                rmatrixgemm<Precision>(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(s2, n, k, alpha, a, ia, ja+s1, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
            }
            return;
        }
        if( n>=m && n>=k )
        {
            
            //
            // A*B = A*(B1 B2)
            //
            ablassplitlength<Precision>(a, n, s1, s2);
            if( optypeb==0 )
            {
                rmatrixgemm<Precision>(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, s2, k, alpha, a, ia, ja, optypea, b, ib, jb+s1, optypeb, beta, c, ic, jc+s1);
            }
            else
            {
                rmatrixgemm<Precision>(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, s2, k, alpha, a, ia, ja, optypea, b, ib+s1, jb, optypeb, beta, c, ic, jc+s1);
            }
            return;
        }
        if( k>=m && k>=n )
        {
            
            //
            // A*B = (A1 A2)*(B1 B2)^T
            //
            ablassplitlength<Precision>(a, k, s1, s2);
            if( optypea==0 && optypeb==0 )
            {
                rmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib+s1, jb, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea==0 && optypeb!=0 )
            {
                rmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib, jb+s1, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea!=0 && optypeb==0 )
            {
                rmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib+s1, jb, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            if( optypea!=0 && optypeb!=0 )
            {
                rmatrixgemm<Precision>(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
                rmatrixgemm<Precision>(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib, jb+s1, optypeb, amp::ampf<Precision>("1.0"), c, ic, jc);
            }
            return;
        }
    }


    /*************************************************************************
    Complex ABLASSplitLength

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ablasinternalsplitlength(int n,
        int nb,
        int& n1,
        int& n2)
    {
        int r;


        if( n<=nb )
        {
            
            //
            // Block size, no further splitting
            //
            n1 = n;
            n2 = 0;
        }
        else
        {
            
            //
            // Greater than block size
            //
            if( n%nb!=0 )
            {
                
                //
                // Split remainder
                //
                n2 = n%nb;
                n1 = n-n2;
            }
            else
            {
                
                //
                // Split on block boundaries
                //
                n2 = n/2;
                n1 = n-n2;
                if( n1%nb==0 )
                {
                    return;
                }
                r = nb-n1%nb;
                n1 = n1+r;
                n2 = n2-r;
            }
        }
    }


    /*************************************************************************
    Level 2 variant of CMatrixRightTRSM
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrighttrsm2(int m,
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
        int i;
        int j;
        amp::campf<Precision> vc;
        amp::campf<Precision> vd;
        int i_;
        int i1_;


        
        //
        // Special case
        //
        if( n*m==0 )
        {
            return;
        }
        
        //
        // Try to call fast TRSM
        //
        if( ablasf::cmatrixrighttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
        {
            return;
        }
        
        //
        // General case
        //
        if( isupper )
        {
            
            //
            // Upper triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // X*A^(-1)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( isunit )
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                        if( j<n-1 )
                        {
                            vc = x(i2+i,j2+j);
                            i1_ = (j1+j+1) - (j2+j+1);
                            for(i_=j2+j+1; i_<=j2+n-1;i_++)
                            {
                                x(i2+i,i_) = x(i2+i,i_) - vc*a(i1+j,i_+i1_);
                            }
                        }
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // X*A^(-T)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=n-1; j>=0; j--)
                    {
                        vc = 0;
                        vd = 1;
                        if( j<n-1 )
                        {
                            i1_ = (j1+j+1)-(j2+j+1);
                            vc = 0.0;
                            for(i_=j2+j+1; i_<=j2+n-1;i_++)
                            {
                                vc += x(i2+i,i_)*a(i1+j,i_+i1_);
                            }
                        }
                        if( !isunit )
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                    }
                }
                return;
            }
            if( optype==2 )
            {
                
                //
                // X*A^(-H)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=n-1; j>=0; j--)
                    {
                        vc = 0;
                        vd = 1;
                        if( j<n-1 )
                        {
                            i1_ = (j1+j+1)-(j2+j+1);
                            vc = 0.0;
                            for(i_=j2+j+1; i_<=j2+n-1;i_++)
                            {
                                vc += x(i2+i,i_)*amp::conj(a(i1+j,i_+i1_));
                            }
                        }
                        if( !isunit )
                        {
                            vd = amp::conj<Precision>(a(i1+j,j1+j));
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                    }
                }
                return;
            }
        }
        else
        {
            
            //
            // Lower triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // X*A^(-1)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=n-1; j>=0; j--)
                    {
                        if( isunit )
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                        if( j>0 )
                        {
                            vc = x(i2+i,j2+j);
                            i1_ = (j1) - (j2);
                            for(i_=j2; i_<=j2+j-1;i_++)
                            {
                                x(i2+i,i_) = x(i2+i,i_) - vc*a(i1+j,i_+i1_);
                            }
                        }
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // X*A^(-T)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        vc = 0;
                        vd = 1;
                        if( j>0 )
                        {
                            i1_ = (j1)-(j2);
                            vc = 0.0;
                            for(i_=j2; i_<=j2+j-1;i_++)
                            {
                                vc += x(i2+i,i_)*a(i1+j,i_+i1_);
                            }
                        }
                        if( !isunit )
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                    }
                }
                return;
            }
            if( optype==2 )
            {
                
                //
                // X*A^(-H)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        vc = 0;
                        vd = 1;
                        if( j>0 )
                        {
                            i1_ = (j1)-(j2);
                            vc = 0.0;
                            for(i_=j2; i_<=j2+j-1;i_++)
                            {
                                vc += x(i2+i,i_)*amp::conj(a(i1+j,i_+i1_));
                            }
                        }
                        if( !isunit )
                        {
                            vd = amp::conj<Precision>(a(i1+j,j1+j));
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                    }
                }
                return;
            }
        }
    }


    /*************************************************************************
    Level-2 subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlefttrsm2(int m,
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
        int i;
        int j;
        amp::campf<Precision> vc;
        amp::campf<Precision> vd;
        int i_;


        
        //
        // Special case
        //
        if( n*m==0 )
        {
            return;
        }
        
        //
        // Try to call fast TRSM
        //
        if( ablasf::cmatrixlefttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
        {
            return;
        }
        
        //
        // General case
        //
        if( isupper )
        {
            
            //
            // Upper triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // A^(-1)*X
                //
                for(i=m-1; i>=0; i--)
                {
                    for(j=i+1; j<=m-1; j++)
                    {
                        vc = a(i1+i,j1+j);
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+i,i_) = x(i2+i,i_) - vc*x(i2+j,i_);
                        }
                    }
                    if( !isunit )
                    {
                        vd = 1/a(i1+i,j1+i);
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+i,i_) = vd*x(i2+i,i_);
                        }
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // A^(-T)*X
                //
                for(i=0; i<=m-1; i++)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+i,j1+i);
                    }
                    for(i_=j2; i_<=j2+n-1;i_++)
                    {
                        x(i2+i,i_) = vd*x(i2+i,i_);
                    }
                    for(j=i+1; j<=m-1; j++)
                    {
                        vc = a(i1+i,j1+j);
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+j,i_) = x(i2+j,i_) - vc*x(i2+i,i_);
                        }
                    }
                }
                return;
            }
            if( optype==2 )
            {
                
                //
                // A^(-H)*X
                //
                for(i=0; i<=m-1; i++)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/amp::conj<Precision>(a(i1+i,j1+i));
                    }
                    for(i_=j2; i_<=j2+n-1;i_++)
                    {
                        x(i2+i,i_) = vd*x(i2+i,i_);
                    }
                    for(j=i+1; j<=m-1; j++)
                    {
                        vc = amp::conj<Precision>(a(i1+i,j1+j));
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+j,i_) = x(i2+j,i_) - vc*x(i2+i,i_);
                        }
                    }
                }
                return;
            }
        }
        else
        {
            
            //
            // Lower triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // A^(-1)*X
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=i-1; j++)
                    {
                        vc = a(i1+i,j1+j);
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+i,i_) = x(i2+i,i_) - vc*x(i2+j,i_);
                        }
                    }
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+j,j1+j);
                    }
                    for(i_=j2; i_<=j2+n-1;i_++)
                    {
                        x(i2+i,i_) = vd*x(i2+i,i_);
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // A^(-T)*X
                //
                for(i=m-1; i>=0; i--)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+i,j1+i);
                    }
                    for(i_=j2; i_<=j2+n-1;i_++)
                    {
                        x(i2+i,i_) = vd*x(i2+i,i_);
                    }
                    for(j=i-1; j>=0; j--)
                    {
                        vc = a(i1+i,j1+j);
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+j,i_) = x(i2+j,i_) - vc*x(i2+i,i_);
                        }
                    }
                }
                return;
            }
            if( optype==2 )
            {
                
                //
                // A^(-H)*X
                //
                for(i=m-1; i>=0; i--)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/amp::conj<Precision>(a(i1+i,j1+i));
                    }
                    for(i_=j2; i_<=j2+n-1;i_++)
                    {
                        x(i2+i,i_) = vd*x(i2+i,i_);
                    }
                    for(j=i-1; j>=0; j--)
                    {
                        vc = amp::conj<Precision>(a(i1+i,j1+j));
                        for(i_=j2; i_<=j2+n-1;i_++)
                        {
                            x(i2+j,i_) = x(i2+j,i_) - vc*x(i2+i,i_);
                        }
                    }
                }
                return;
            }
        }
    }


    /*************************************************************************
    Level 2 subroutine

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrighttrsm2(int m,
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
        int i;
        int j;
        amp::ampf<Precision> vr;
        amp::ampf<Precision> vd;


        
        //
        // Special case
        //
        if( n*m==0 )
        {
            return;
        }
        
        //
        // Try to use "fast" code
        //
        if( ablasf::rmatrixrighttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
        {
            return;
        }
        
        //
        // General case
        //
        if( isupper )
        {
            
            //
            // Upper triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // X*A^(-1)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( isunit )
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                        if( j<n-1 )
                        {
                            vr = x(i2+i,j2+j);
                            amp::vsub(x.getrow(i2+i, j2+j+1, j2+n-1), a.getrow(i1+j, j1+j+1, j1+n-1), vr);
                        }
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // X*A^(-T)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=n-1; j>=0; j--)
                    {
                        vr = 0;
                        vd = 1;
                        if( j<n-1 )
                        {
                            vr = amp::vdotproduct(x.getrow(i2+i, j2+j+1, j2+n-1), a.getrow(i1+j, j1+j+1, j1+n-1));
                        }
                        if( !isunit )
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vr)/vd;
                    }
                }
                return;
            }
        }
        else
        {
            
            //
            // Lower triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // X*A^(-1)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=n-1; j>=0; j--)
                    {
                        if( isunit )
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                        if( j>0 )
                        {
                            vr = x(i2+i,j2+j);
                            amp::vsub(x.getrow(i2+i, j2, j2+j-1), a.getrow(i1+j, j1, j1+j-1), vr);
                        }
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // X*A^(-T)
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        vr = 0;
                        vd = 1;
                        if( j>0 )
                        {
                            vr = amp::vdotproduct(x.getrow(i2+i, j2, j2+j-1), a.getrow(i1+j, j1, j1+j-1));
                        }
                        if( !isunit )
                        {
                            vd = a(i1+j,j1+j);
                        }
                        x(i2+i,j2+j) = (x(i2+i,j2+j)-vr)/vd;
                    }
                }
                return;
            }
        }
    }


    /*************************************************************************
    Level 2 subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlefttrsm2(int m,
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
        int i;
        int j;
        amp::ampf<Precision> vr;
        amp::ampf<Precision> vd;


        
        //
        // Special case
        //
        if( n*m==0 )
        {
            return;
        }
        
        //
        // Try fast code
        //
        if( ablasf::rmatrixlefttrsmf<Precision>(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
        {
            return;
        }
        
        //
        // General case
        //
        if( isupper )
        {
            
            //
            // Upper triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // A^(-1)*X
                //
                for(i=m-1; i>=0; i--)
                {
                    for(j=i+1; j<=m-1; j++)
                    {
                        vr = a(i1+i,j1+j);
                        amp::vsub(x.getrow(i2+i, j2, j2+n-1), x.getrow(i2+j, j2, j2+n-1), vr);
                    }
                    if( !isunit )
                    {
                        vd = 1/a(i1+i,j1+i);
                        amp::vmul(x.getrow(i2+i, j2, j2+n-1), vd);
                    }
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // A^(-T)*X
                //
                for(i=0; i<=m-1; i++)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+i,j1+i);
                    }
                    amp::vmul(x.getrow(i2+i, j2, j2+n-1), vd);
                    for(j=i+1; j<=m-1; j++)
                    {
                        vr = a(i1+i,j1+j);
                        amp::vsub(x.getrow(i2+j, j2, j2+n-1), x.getrow(i2+i, j2, j2+n-1), vr);
                    }
                }
                return;
            }
        }
        else
        {
            
            //
            // Lower triangular matrix
            //
            if( optype==0 )
            {
                
                //
                // A^(-1)*X
                //
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=i-1; j++)
                    {
                        vr = a(i1+i,j1+j);
                        amp::vsub(x.getrow(i2+i, j2, j2+n-1), x.getrow(i2+j, j2, j2+n-1), vr);
                    }
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+j,j1+j);
                    }
                    amp::vmul(x.getrow(i2+i, j2, j2+n-1), vd);
                }
                return;
            }
            if( optype==1 )
            {
                
                //
                // A^(-T)*X
                //
                for(i=m-1; i>=0; i--)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = 1/a(i1+i,j1+i);
                    }
                    amp::vmul(x.getrow(i2+i, j2, j2+n-1), vd);
                    for(j=i-1; j>=0; j--)
                    {
                        vr = a(i1+i,j1+j);
                        amp::vsub(x.getrow(i2+j, j2, j2+n-1), x.getrow(i2+i, j2, j2+n-1), vr);
                    }
                }
                return;
            }
        }
    }


    /*************************************************************************
    Level 2 subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixsyrk2(int n,
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
        int i;
        int j;
        int j1;
        int j2;
        amp::campf<Precision> v;
        int i_;
        int i1_;


        
        //
        // Fast exit (nothing to be done)
        //
        if( (alpha==0 || k==0) && beta==1 )
        {
            return;
        }
        
        //
        // Try to call fast SYRK
        //
        if( ablasf::cmatrixsyrkf<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper) )
        {
            return;
        }
        
        //
        // SYRK
        //
        if( optypea==0 )
        {
            
            //
            // C=alpha*A*A^H+beta*C
            //
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = i;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = i;
                }
                for(j=j1; j<=j2; j++)
                {
                    if( alpha!=0 && k>0 )
                    {
                        v = 0.0;
                        for(i_=ja; i_<=ja+k-1;i_++)
                        {
                            v += a(ia+i,i_)*amp::conj(a(ia+j,i_));
                        }
                    }
                    else
                    {
                        v = 0;
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        else
        {
            
            //
            // C=alpha*A^H*A+beta*C
            //
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = i;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = i;
                }
                if( beta==0 )
                {
                    for(j=j1; j<=j2; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
                else
                {
                    for(i_=jc+j1; i_<=jc+j2;i_++)
                    {
                        c(ic+i,i_) = beta*c(ic+i,i_);
                    }
                }
            }
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( isupper )
                    {
                        j1 = j;
                        j2 = n-1;
                    }
                    else
                    {
                        j1 = 0;
                        j2 = j;
                    }
                    v = alpha*amp::conj<Precision>(a(ia+i,ja+j));
                    i1_ = (ja+j1) - (jc+j1);
                    for(i_=jc+j1; i_<=jc+j2;i_++)
                    {
                        c(ic+j,i_) = c(ic+j,i_) + v*a(ia+i,i_+i1_);
                    }
                }
            }
            return;
        }
    }


    /*************************************************************************
    Level 2 subrotuine
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsyrk2(int n,
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
        int i;
        int j;
        int j1;
        int j2;
        amp::ampf<Precision> v;


        
        //
        // Fast exit (nothing to be done)
        //
        if( (alpha==0 || k==0) && beta==1 )
        {
            return;
        }
        
        //
        // Try to call fast SYRK
        //
        if( ablasf::rmatrixsyrkf<Precision>(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper) )
        {
            return;
        }
        
        //
        // SYRK
        //
        if( optypea==0 )
        {
            
            //
            // C=alpha*A*A^H+beta*C
            //
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = i;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = i;
                }
                for(j=j1; j<=j2; j++)
                {
                    if( alpha!=0 && k>0 )
                    {
                        v = amp::vdotproduct(a.getrow(ia+i, ja, ja+k-1), a.getrow(ia+j, ja, ja+k-1));
                    }
                    else
                    {
                        v = 0;
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        else
        {
            
            //
            // C=alpha*A^H*A+beta*C
            //
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = i;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = i;
                }
                if( beta==0 )
                {
                    for(j=j1; j<=j2; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
                else
                {
                    amp::vmul(c.getrow(ic+i, jc+j1, jc+j2), beta);
                }
            }
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( isupper )
                    {
                        j1 = j;
                        j2 = n-1;
                    }
                    else
                    {
                        j1 = 0;
                        j2 = j;
                    }
                    v = alpha*a(ia+i,ja+j);
                    amp::vadd(c.getrow(ic+j, jc+j1, jc+j2), a.getrow(ia+i, ja+j1, ja+j2), v);
                }
            }
            return;
        }
    }


    /*************************************************************************
    GEMM kernel

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixgemmk(int m,
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
        int i;
        int j;
        amp::campf<Precision> v;
        int i_;
        int i1_;


        
        //
        // Special case
        //
        if( m*n==0 )
        {
            return;
        }
        
        //
        // Try optimized code
        //
        if( ablasf::cmatrixgemmf<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc) )
        {
            return;
        }
        
        //
        // Another special case
        //
        if( k==0 )
        {
            if( beta!=0 )
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j);
                    }
                }
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
            }
            return;
        }
        
        //
        // General case
        //
        if( optypea==0 && optypeb!=0 )
        {
            
            //
            // A*B'
            //
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( k==0 || alpha==0 )
                    {
                        v = 0;
                    }
                    else
                    {
                        if( optypeb==1 )
                        {
                            i1_ = (jb)-(ja);
                            v = 0.0;
                            for(i_=ja; i_<=ja+k-1;i_++)
                            {
                                v += a(ia+i,i_)*b(ib+j,i_+i1_);
                            }
                        }
                        else
                        {
                            i1_ = (jb)-(ja);
                            v = 0.0;
                            for(i_=ja; i_<=ja+k-1;i_++)
                            {
                                v += a(ia+i,i_)*amp::conj(b(ib+j,i_+i1_));
                            }
                        }
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        if( optypea==0 && optypeb==0 )
        {
            
            //
            // A*B
            //
            for(i=0; i<=m-1; i++)
            {
                if( beta!=0 )
                {
                    for(i_=jc; i_<=jc+n-1;i_++)
                    {
                        c(ic+i,i_) = beta*c(ic+i,i_);
                    }
                }
                else
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
                if( alpha!=0 )
                {
                    for(j=0; j<=k-1; j++)
                    {
                        v = alpha*a(ia+i,ja+j);
                        i1_ = (jb) - (jc);
                        for(i_=jc; i_<=jc+n-1;i_++)
                        {
                            c(ic+i,i_) = c(ic+i,i_) + v*b(ib+j,i_+i1_);
                        }
                    }
                }
            }
            return;
        }
        if( optypea!=0 && optypeb!=0 )
        {
            
            //
            // A'*B'
            //
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( alpha==0 )
                    {
                        v = 0;
                    }
                    else
                    {
                        if( optypea==1 )
                        {
                            if( optypeb==1 )
                            {
                                i1_ = (jb)-(ia);
                                v = 0.0;
                                for(i_=ia; i_<=ia+k-1;i_++)
                                {
                                    v += a(i_,ja+i)*b(ib+j,i_+i1_);
                                }
                            }
                            else
                            {
                                i1_ = (jb)-(ia);
                                v = 0.0;
                                for(i_=ia; i_<=ia+k-1;i_++)
                                {
                                    v += a(i_,ja+i)*amp::conj(b(ib+j,i_+i1_));
                                }
                            }
                        }
                        else
                        {
                            if( optypeb==1 )
                            {
                                i1_ = (jb)-(ia);
                                v = 0.0;
                                for(i_=ia; i_<=ia+k-1;i_++)
                                {
                                    v += amp::conj(a(i_,ja+i))*b(ib+j,i_+i1_);
                                }
                            }
                            else
                            {
                                i1_ = (jb)-(ia);
                                v = 0.0;
                                for(i_=ia; i_<=ia+k-1;i_++)
                                {
                                    v += amp::conj(a(i_,ja+i))*amp::conj(b(ib+j,i_+i1_));
                                }
                            }
                        }
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        if( optypea!=0 && optypeb==0 )
        {
            
            //
            // A'*B
            //
            if( beta==0 )
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    for(i_=jc; i_<=jc+n-1;i_++)
                    {
                        c(ic+i,i_) = beta*c(ic+i,i_);
                    }
                }
            }
            if( alpha!=0 )
            {
                for(j=0; j<=k-1; j++)
                {
                    for(i=0; i<=m-1; i++)
                    {
                        if( optypea==1 )
                        {
                            v = alpha*a(ia+j,ja+i);
                        }
                        else
                        {
                            v = alpha*amp::conj<Precision>(a(ia+j,ja+i));
                        }
                        i1_ = (jb) - (jc);
                        for(i_=jc; i_<=jc+n-1;i_++)
                        {
                            c(ic+i,i_) = c(ic+i,i_) + v*b(ib+j,i_+i1_);
                        }
                    }
                }
            }
            return;
        }
    }


    /*************************************************************************
    GEMM kernel

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixgemmk(int m,
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
        int i;
        int j;
        amp::ampf<Precision> v;


        
        //
        // if matrix size is zero
        //
        if( m*n==0 )
        {
            return;
        }
        
        //
        // Try optimized code
        //
        if( ablasf::rmatrixgemmf<Precision>(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc) )
        {
            return;
        }
        
        //
        // if K=0, then C=Beta*C
        //
        if( k==0 )
        {
            if( beta!=1 )
            {
                if( beta!=0 )
                {
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            c(ic+i,jc+j) = beta*c(ic+i,jc+j);
                        }
                    }
                }
                else
                {
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            c(ic+i,jc+j) = 0;
                        }
                    }
                }
            }
            return;
        }
        
        //
        // General case
        //
        if( optypea==0 && optypeb!=0 )
        {
            
            //
            // A*B'
            //
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( k==0 || alpha==0 )
                    {
                        v = 0;
                    }
                    else
                    {
                        v = amp::vdotproduct(a.getrow(ia+i, ja, ja+k-1), b.getrow(ib+j, jb, jb+k-1));
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        if( optypea==0 && optypeb==0 )
        {
            
            //
            // A*B
            //
            for(i=0; i<=m-1; i++)
            {
                if( beta!=0 )
                {
                    amp::vmul(c.getrow(ic+i, jc, jc+n-1), beta);
                }
                else
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
                if( alpha!=0 )
                {
                    for(j=0; j<=k-1; j++)
                    {
                        v = alpha*a(ia+i,ja+j);
                        amp::vadd(c.getrow(ic+i, jc, jc+n-1), b.getrow(ib+j, jb, jb+n-1), v);
                    }
                }
            }
            return;
        }
        if( optypea!=0 && optypeb!=0 )
        {
            
            //
            // A'*B'
            //
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( alpha==0 )
                    {
                        v = 0;
                    }
                    else
                    {
                        v = amp::vdotproduct(a.getcolumn(ja+i, ia, ia+k-1), b.getrow(ib+j, jb, jb+k-1));
                    }
                    if( beta==0 )
                    {
                        c(ic+i,jc+j) = alpha*v;
                    }
                    else
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                    }
                }
            }
            return;
        }
        if( optypea!=0 && optypeb==0 )
        {
            
            //
            // A'*B
            //
            if( beta==0 )
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    amp::vmul(c.getrow(ic+i, jc, jc+n-1), beta);
                }
            }
            if( alpha!=0 )
            {
                for(j=0; j<=k-1; j++)
                {
                    for(i=0; i<=m-1; i++)
                    {
                        v = alpha*a(ia+j,ja+i);
                        amp::vadd(c.getrow(ic+i, jc, jc+n-1), b.getrow(ib+j, jb, jb+n-1), v);
                    }
                }
            }
            return;
        }
    }
} // namespace

#endif
