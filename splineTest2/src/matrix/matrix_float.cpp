/*=============================================================================
        File: matrix.cpp
     Purpose:       
    Revision: $Id: matrix_float.cpp,v 1.2 2002/05/13 21:07:45 philosophil Exp $
  Created by: Philippe Lavoie          (3 Oct, 1996)
 Modified by: 

 Copyright notice:
          Copyright (C) 1996-1998 Philippe Lavoie
 
	  This library is free software; you can redistribute it and/or
	  modify it under the terms of the GNU Library General Public
	  License as published by the Free Software Foundation; either
	  version 2 of the License, or (at your option) any later version.
 
	  This library is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  Library General Public License for more details.
 
	  You should have received a copy of the GNU Library General Public
	  License along with this library; if not, write to the Free
	  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
=============================================================================*/

#include "matrix.cpp"

namespace PLib {
  
  template<> void Matrix<float>::qSort(){
    qsort((char*)m,rows()*cols(),sizeof(float),compareFloat) ;
  }

#ifdef NO_IMPLICIT_TEMPLATES

  template class Matrix<float> ;
  
  template Matrix<float> operator+(const Matrix<float>&,const Matrix<float>&);
  template Matrix<float> operator-(const Matrix<float>&,const Matrix<float>&);
  template Matrix<float> operator*(const Matrix<float>&,const Matrix<float>&);
  template Matrix<float> operator*(const double,const Matrix<float>&);
  template Matrix<float> operator*(const Complex&,const Matrix<float>&);
  template Vector<float> operator*(const Matrix<float>&,const Vector<float>&);
  template int operator==(const Matrix<float>&,const Matrix<float>&);
  //template int operator!=(const Matrix<float>&,const Matrix<float>&);
  template Matrix<float> comm(const Matrix<float>&,const Matrix<float>&);


#endif

}
