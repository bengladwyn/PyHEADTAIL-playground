/*************************************************************************
AMP library for multiple precision operations.
Copyright (c) 2007-2009 Sergey Bochkanov (ALGLIB project).

>>> LICENSE >>>
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
#include "stdafx.h"
#include "amp.h"

/************************************************************************
Storage of mpfr lists
************************************************************************/
std::vector< amp::mpfr_record* > _mpfr_storage_v;
gmp_randstate_t _mpfr_storage_rs;
bool _mpfr_storage_init = false;

amp::mpfr_record* amp::mpfr_storage::newMpfr(unsigned int Precision)
{
    amp::mpfr_record_ptr &lst = getList(Precision);
    if( lst==NULL )
    {
        amp::mpfr_record *rec = new amp::mpfr_record;
        rec->Precision = Precision;
        rec->refCount = 0;
        mpfr_init2(rec->value, Precision);
        rec->next = lst;
        lst = rec;
    }
    
    amp::mpfr_record *p = lst;
    p->refCount = 1;
    lst = lst->next;
    return p;
}

void amp::mpfr_storage::deleteMpfr(amp::mpfr_record* ref)
{
    amp::mpfr_record_ptr &lst = getList(ref->Precision);
    ref->next = lst;
    lst = ref;
}

gmp_randstate_t* amp::mpfr_storage::getRandState()
{
    if( !_mpfr_storage_init )
    {
        time_t _timer;
        gmp_randinit_default(_mpfr_storage_rs);
        gmp_randseed_ui(_mpfr_storage_rs, (unsigned long int)(time(&_timer)));
        _mpfr_storage_init = true;
    }
    return &_mpfr_storage_rs;
}

void amp::mpfr_storage::seedRandState(unsigned long seed)
{
    //
    // call getRandState to unsure that state is initialized
    //
    getRandState();

    //
    // seed it
    //
    gmp_randseed_ui(_mpfr_storage_rs, (unsigned long int)seed);
}

/*amp::mpfr_storage::clearStorage()
{
    unsigned int i;
    amp::mpfr_record *p;
    for(i=0; i<v.size(); i++)
        while( v[i]!=NULL )
        {
            p = _mpfr_storage_v[i]->next;
            mpfr_clear(_mpfr_storage_v[i]->value);
            delete _mpfr_storage_v[i];
            _mpfr_storage_v[i] = p;
        }
}*/

amp::mpfr_record_ptr& amp::mpfr_storage::getList(unsigned int Precision)
{
    static amp::mpfr_record_ptr tmp       = NULL;
    static unsigned int lastPrec          = -1;
    static amp::mpfr_record_ptr &lastList = tmp;
    if( lastPrec!=Precision )
    {
        while( _mpfr_storage_v.size()<Precision+1 )
            _mpfr_storage_v.push_back(NULL);
        lastPrec = Precision;
        lastList = _mpfr_storage_v[Precision];
    }
    return lastList;
}


/************************************************************************
Storage of mpfr lists
************************************************************************/
amp::mpfr_reference::mpfr_reference()
{
    ref = NULL;
}

amp::mpfr_reference::mpfr_reference(const amp::mpfr_reference& r)
{
    ref = r.ref;
    if( ref!=NULL )
        ref->refCount++;
}

amp::mpfr_reference& amp::mpfr_reference::operator= (
    const amp::mpfr_reference &r)
{
    if ( &r==this )
        return *this;
    if ( ref==r.ref )
        return *this;
    if( ref!=NULL )
        free();
    ref = r.ref;
    if( ref!=NULL )
        ref->refCount++;
    return *this;
}

amp::mpfr_reference::~mpfr_reference()
{
    if( ref!=NULL )
        free();
}
        
void amp::mpfr_reference::initialize(int Precision)
{
    if( ref!=NULL )
        free();
    ref = amp::mpfr_storage::newMpfr(Precision);
    ref->refCount = 1;
}

void amp::mpfr_reference::free()
{
    if( ref==NULL )
        throw amp::internalError();
    ref->refCount--;
    if( ref->refCount==0 )
        amp::mpfr_storage::deleteMpfr(ref);
    ref = NULL;
}
        
mpfr_srcptr amp::mpfr_reference::getReadPtr() const
{
    if( ref==NULL )
        throw amp::internalError();
    return ref->value;
}

mpfr_ptr amp::mpfr_reference::getWritePtr()
{
    if( ref==NULL )
        throw amp::internalError();
    if( ref->refCount==1 )
        return ref->value;

    amp::mpfr_record *newref = amp::mpfr_storage::newMpfr(ref->Precision);
    mpfr_set(newref->value, ref->value, GMP_RNDN);
    
    free();
    ref = newref;
    return ref->value;
}
