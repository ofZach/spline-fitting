#include "testApp.h"




#include <string.h>
#include <assert.h>


const double EPSILON_VALUE = 4.37114e-05;



template<typename T>
BandedMatrix<T>::BandedMatrix( int iSize, int iLBands, int iUBands )
{
	assert(iSize > 0 && iLBands >= 0 && iUBands >= 0);
	assert(iLBands < iSize && iUBands < iSize);
    
	m_iSize = iSize;
	m_iLBands = iLBands;
	m_iUBands = iUBands;
	allocate();
}

template<typename T>
BandedMatrix<T>::BandedMatrix( const BandedMatrix& rkM )
{
	m_afDBand = 0;
	m_aafLBand = 0;
	m_aafUBand = 0;
	*this = rkM;
}

template<typename T>
BandedMatrix<T>::~BandedMatrix()
{
	deallocate();
}

template<typename T>
BandedMatrix<T>& BandedMatrix<T>::operator=( const BandedMatrix& rkM)
{
	deallocate();
	m_iSize = rkM.m_iSize;
	m_iLBands = rkM.m_iLBands;
	m_iUBands = rkM.m_iUBands;
	allocate();
    
	size_t uiSize = m_iSize*sizeof(T);
	memcpy( m_afDBand, rkM.m_afDBand, uiSize );
    
	int i;
	for (i = 0; i < m_iLBands; i++) {
		uiSize = (m_iSize-1-i)*sizeof(T);
		memcpy( m_aafLBand[i], rkM.m_aafLBand[i], uiSize );
	}
    
	for (i = 0; i < m_iUBands; i++) {
		uiSize = (m_iSize-1-i)*sizeof(T);
		memcpy( m_aafUBand[i], rkM.m_aafUBand[i], uiSize );
	}
    
	return *this;
}

template<typename T>
int BandedMatrix<T>::getSize() const
{
    return m_iSize;
}

template<typename T>
int BandedMatrix<T>::getLBands() const
{
    return m_iLBands;
}

template<typename T>
int BandedMatrix<T>::getUBands() const
{
    return m_iUBands;
}

template<typename T>
T* BandedMatrix<T>::getDBand()
{
    return m_afDBand;
}

template<typename T>
const T* BandedMatrix<T>::getDBand() const
{
    return m_afDBand;
}

template<typename T>
int BandedMatrix<T>::getLBandMax( int i ) const
{
    assert(0 <= i && i < m_iLBands);
    return m_iSize-1-i;
}

template<typename T>
T* BandedMatrix<T>::getLBand( int i )
{
    if ( m_aafLBand )
    {
        assert(0 <= i && i < m_iLBands);
        return m_aafLBand[i];
    }
    return 0;
}

template<typename T>
const T* BandedMatrix<T>::getLBand( int i ) const
{
    if (m_aafLBand)
    {
        assert(0 <= i && i < m_iLBands);
        return m_aafLBand[i];
    }
    return 0;
}

template<typename T>
int BandedMatrix<T>::getUBandMax( int i ) const
{
    assert(0 <= i && i < m_iUBands);
    return m_iSize-1-i;
}

template<typename T>
T* BandedMatrix<T>::getUBand( int i )
{
    if (m_aafUBand)
    {
        assert(0 <= i && i < m_iUBands);
        return m_aafUBand[i];
    }
    return 0;
}

template<typename T>
const T* BandedMatrix<T>::getUBand( int i ) const
{
    if (m_aafUBand)
    {
        assert(0 <= i && i < m_iUBands);
        return m_aafUBand[i];
    }
    return 0;
}

template<typename T>
T& BandedMatrix<T>::operator()( int iRow, int iCol )
{
    assert(0 <= iRow && iRow < m_iSize && 0 <= iCol && iCol < m_iSize);
    
    int iBand = iCol - iRow;
    if (iBand > 0)
    {
        if (--iBand < m_iUBands && iRow < m_iSize-1-iBand)
        {
            return m_aafUBand[iBand][iRow];
        }
    }
    else if ( iBand < 0 )
    {
        iBand = -iBand;
        if (--iBand < m_iLBands && iCol < m_iSize-1-iBand)
        {
            return m_aafLBand[iBand][iCol];
        }
    }
    else
    {
        return m_afDBand[iRow];
    }
    
    static T s_fDummy = (T)0.0;
    return s_fDummy;
}

template<typename T>
T BandedMatrix<T>::operator()( int iRow, int iCol ) const
{
    assert(0 <= iRow && iRow < m_iSize && 0 <= iCol && iCol < m_iSize);
    
    int iBand = iCol - iRow;
    if (iBand > 0)
    {
        if (--iBand < m_iUBands && iRow < m_iSize-1-iBand)
        {
            return m_aafUBand[iBand][iRow];
        }
    }
    else if ( iBand < 0 )
    {
        iBand = -iBand;
        if (--iBand < m_iLBands && iCol < m_iSize-1-iBand)
        {
            return m_aafLBand[iBand][iCol];
        }
    }
    else
    {
        return m_afDBand[iRow];
    }
    
    return (T)0.0;
}

template<typename T>
void BandedMatrix<T>::setZero()
{
    assert(m_iSize > 0);
    
    memset(m_afDBand,0,m_iSize*sizeof(T));
    
    int i;
    for (i = 0; i < m_iLBands; i++)
    {
        memset(m_aafLBand[i],0,(m_iSize-1-i)*sizeof(T));
    }
    
    for (i = 0; i < m_iUBands; i++)
    {
        memset(m_aafUBand[i],0,(m_iSize-1-i)*sizeof(T));
    }
}

template<typename T>
void BandedMatrix<T>::setIdentity()
{
    assert(m_iSize > 0);
    
    int i;
    for (i = 0; i < m_iSize; i++)
    {
        m_afDBand[i] = (T)1.0;
    }
    
    for (i = 0; i < m_iLBands; i++)
    {
        memset(m_aafLBand[i],0,(m_iSize-1-i)*sizeof(T));
    }
    
    for (i = 0; i < m_iUBands; i++)
    {
        memset(m_aafUBand[i],0,(m_iSize-1-i)*sizeof(T));
    }
}

template<typename T>
void BandedMatrix<T>::allocate()
{
    // assert:  m_iSize, m_iLBands, m_iRBandQuantity already set
    // assert:  m_afDBand, m_aafLBand, m_aafUBand all null
    
    m_afDBand = new T[m_iSize];
    memset( m_afDBand,0,m_iSize*sizeof(T) );
    
    if( m_iLBands > 0 )
    {
        m_aafLBand = new T*[m_iLBands];
    }
    else
    {
        m_aafLBand = 0;
    }
    
    if (m_iUBands > 0)
    {
        m_aafUBand = new T*[m_iUBands];
    }
    else
    {
        m_aafUBand = 0;
    }
    
    int i;
    for( i = 0; i < m_iLBands; i++ ) {
        m_aafLBand[i] = new T[m_iSize-1-i];
        memset(m_aafLBand[i],0,(m_iSize-1-i)*sizeof(T));
    }
    
    for (i = 0; i < m_iUBands; i++) {
        m_aafUBand[i] = new T[m_iSize-1-i];
        memset( m_aafUBand[i],0,(m_iSize-1-i)*sizeof(T) );
    }
}

template<typename T>
void BandedMatrix<T>::deallocate()
{
    delete [] m_afDBand;
    
    int i;
    
    if( m_aafLBand ) {
        for (i = 0; i < m_iLBands; i++)
        {
            delete [] m_aafLBand[i];
        }
        
        delete [] m_aafLBand;
        m_aafLBand = 0;
    }
    
    if( m_aafUBand ) {
        for (i = 0; i < m_iUBands; i++) {
            delete [] m_aafUBand[i];
        }
        
        delete [] m_aafUBand;
        m_aafUBand = 0;
    }
}

template class BandedMatrix<float>;
template class BandedMatrix<double>;


class BSplineBasis
{
public:
	BSplineBasis();
    
	// Open uniform or periodic uniform.  The knot array is internally
	// generated with equally spaced elements.
	BSplineBasis( int aNumCtrlPoints, int iDegree, bool bOpen );
	void create( int aNumCtrlPoints, int iDegree, bool bOpen );
    
	// Open nonuniform.  The knot array must have n-d elements.  The elements
	// must be nondecreasing.  Each element must be in [0,1].  The caller is
	// responsible for deleting afKnot.  An internal copy is made, so to
	// dynamically change knots you must use the setKnot function.
	BSplineBasis( int aNumCtrlPoints, int iDegree, const float* afKnot );
	void create( int aNumCtrlPoints, int iDegree, const float* afKnot );
    
	BSplineBasis( const BSplineBasis &basis );
	BSplineBasis& operator=( const BSplineBasis &basis );
    
	~BSplineBasis();
    
	int getNumControlPoints() const;
	int getDegree() const;
	bool isOpen() const;
	bool isUniform() const;
    
	// The knot values can be changed only if the basis function is nonuniform
	// and the input index is valid (0 <= i <= n-d-1).  If these conditions
	// are not satisfied, getKnot returns MAX_REAL.
	void setKnot( int i, float fKnot );
	float getKnot( int i ) const;
    
	// access basis functions and their derivatives
	float getD0( int i ) const;
	float getD1( int i ) const;
	float getD2( int i ) const;
	float getD3( int i ) const;
    
	// evaluate basis functions and their derivatives
	void compute( float fTime, unsigned int uiOrder, int &riMinIndex, int &riMaxIndex ) const;
    
protected:
	int initialize( int iNumCtrlPoints, int iDegree, bool bOpen );
	float** allocate() const;
	void deallocate( float** aafArray );
    
	// Determine knot index i for which knot[i] <= rfTime < knot[i+1].
	int getKey( float& rfTime ) const;
    
	int mNumCtrlPoints;    // n+1
	int mDegree;           // d
	float *mKnots;          // knot[n+d+2]
	bool mOpen, mUniform;
    
	// Storage for the basis functions and their derivatives first three
	// derivatives.  The basis array is always allocated by the constructor
	// calls.  A derivative basis array is allocated on the first call to a
	// derivative member function.
	float **m_aafBD0;             // bd0[d+1][n+d+1]
	mutable float **m_aafBD1;     // bd1[d+1][n+d+1]
	mutable float **m_aafBD2;     // bd2[d+1][n+d+1]
	mutable float **m_aafBD3;     // bd3[d+1][n+d+1]
};

template<typename T>
class BSpline
{
public:
	// Construction and destruction.  The caller is responsible for deleting
	// the input arrays if they were dynamically allocated.  Internal copies
	// of the arrays are made, so to dynamically change control points or
	// knots you must use the 'setControlPoint', 'getControlPoint', and
	// 'Knot' member functions.
    
	// Uniform spline.  The number of control points is n+1 >= 2.  The degree
	// of the B-spline is d and must satisfy 1 <= d <= n.  The knots are
	// implicitly calculated in [0,1].  If bOpen is 'true', the spline is
	// open and the knots are
	//   t[i] = 0,               0 <= i <= d
	//          (i-d)/(n+1-d),   d+1 <= i <= n
	//          1,               n+1 <= i <= n+d+1
	// If bOpen is 'false', the spline is periodic and the knots are
	//   t[i] = (i-d)/(n+1-d),   0 <= i <= n+d+1
	// If bLoop is 'true', extra control points are added to generate a closed
	// curve.  For an open spline, the control point array is reallocated and
	// one extra control point is added, set to the first control point
	// C[n+1] = C[0].  For a periodic spline, the control point array is
	// reallocated and the first d points are replicated.  In either case the
	// knot array is calculated accordingly.
	BSpline( const std::vector<T> &points, int degree, bool loop, bool open );
    
	// Open, nonuniform spline.  The knot array must have n-d elements.  The
	// elements must be nondecreasing.  Each element must be in [0,1].
	BSpline() : mCtrlPoints( 0 ), mNumCtrlPoints( -1 ) {}
	BSpline( int numControlPoints, const T *controlPoints, int degree, bool loop, const float *knots );
	BSpline( const BSpline &bspline );
	BSpline& operator=( const BSpline &bspline );
    
	~BSpline(); // good bye
    
	int getNumControlPoints() const { return mNumCtrlPoints; }
	int getDegree() const { return mBasis.getDegree(); }
	int getNumSpans() const { return mNumCtrlPoints - mBasis.getDegree(); }
	bool isOpen() const { return mBasis.isOpen(); }
	bool isUniform() const { return mBasis.isUniform(); }
	bool isLoop() const { return mLoop; }
    
	// Control points may be changed at any time.  The input index should be
	// valid (0 <= i <= n).  If it is invalid, getControlPoint returns a
	// vector whose components are all MAX_REAL.
	void setControlPoint( int i, const T &rkCtrl );
	T getControlPoint( int i ) const;
    
	// The knot values can be changed only if the basis function is nonuniform
	// and the input index is valid (0 <= i <= n-d-1).  If these conditions
	// are not satisfied, getKnot returns MAX_REAL.
	void setKnot( int i, float fKnot );
	float getKnot( int i ) const;
    
	// The spline is defined for 0 <= t <= 1.  If a t-value is outside [0,1],
	// an open spline clamps t to [0,1].  That is, if t > 1, t is set to 1;
	// if t < 0, t is set to 0.  A periodic spline wraps to to [0,1].  That
	// is, if t is outside [0,1], then t is set to t-floor(t).
	T getPosition( float t ) const;
	T getDerivative( float t ) const;
	T getSecondDerivative( float t ) const;
	T getThirdDerivative( float t ) const;
	
    //typename T::TYPE getSpeed( float t ) const { return getDerivative( t ).length(); }
    
	float getLength( float fT0, float fT1 ) const;
    
	// If you need position and derivatives at the same time, it is more
	// efficient to call these functions.  Pass the addresses of those
	// quantities whose values you want.  You may pass 0 in any argument
	// whose value you do not want.
	void get( float t, T *position, T *firstDerivative = NULL, T *secondDerivative = NULL, T *thirdDerivative = NULL ) const;
	//! Returns the time associated with an arc length in the range [0,getLength(0,1)]
	float getTime( float length ) const;
    
	// Access the basis function to compute it without control points.  This
	// is useful for least squares fitting of curves.
	BSplineBasis& getBasis();
    
protected:
    // Replicate the necessary number of control points when the create
    // function has bLoop equal to true, in which case the spline curve must
    // be a closed curve.
    void createControl( const T *akCtrlPoint );
    
    int mNumCtrlPoints;
    T *mCtrlPoints;  // ctrl[n+1]
    bool mLoop;
    BSplineBasis mBasis;
    int mReplicate;  // the number of replicated control points
};

typedef BSpline<ofVec2f> BSpline2f;
typedef BSpline<ofVec3f> BSpline3f;





    template<typename T>
    class BSplineFitBasis {
    public:
        // Construction and destruction.  This class is only for open uniform
        // B-spline basis functions.  The input is the number of control points
        // for a B-spline curve using this basis and the degree of that curve.
        BSplineFitBasis( int iQuantity, int iDegree );
        ~BSplineFitBasis();
        
        // Data member access.
        int getQuantity() const;
        int getDegree() const;
        
        // Evaluate the basis functions.  This function fills in the values
        // returned by GetValue(i) for 0 <= i <= degree.  The return indices iMin
        // and iMax are relative to the array of control points.  The GetValue(i)
        // are the coefficients for the control points ctrl[iMin] throught
        // ctrl[iMax] in the curve evaluation (i.e. the curve has local control).
        void compute( T fT, int &iMin, int &iMax ) const;
        T getValue( int i ) const;
        
    private:
        // The number of control points and degree for the curve.
        int m_iQuantity, m_iDegree;
        
        // The storage for knots and basis evaluation.
        mutable T* m_afValue;  // m_afValue[0..degree]
        mutable T* m_afKnot;   // m_afKnot[2*degree]
    };


//////////////////////////////////////////////////////////////////////////////////////////////
// BSplineBasis
template <class T>
void allocate2D( int iCols, int iRows, T**& raatArray )
{
    raatArray = new T*[iRows];
    raatArray[0] = new T[iRows*iCols];
    for( int iRow = 1; iRow < iRows; iRow++ ) {
        raatArray[iRow] = &raatArray[0][iCols*iRow];
    }
}

template <class T>
void deallocate2D( T**& raatArray )
{
    if( raatArray ) {
        delete raatArray[0];
        delete raatArray;
        raatArray = 0;
    }
}

BSplineBasis::BSplineBasis()
: mKnots( 0 ), m_aafBD0( 0 ), m_aafBD1( 0 ), m_aafBD2( 0 ), m_aafBD3( 0 ), mNumCtrlPoints( -1 )
{
}

BSplineBasis::BSplineBasis( int iNumCtrlPoints, int iDegree, bool bOpen )
{
    create( iNumCtrlPoints, iDegree, bOpen );
}

BSplineBasis::BSplineBasis( const BSplineBasis &basis )
: mKnots( 0 ), m_aafBD0( 0 ), m_aafBD1( 0 ), m_aafBD2( 0 ), m_aafBD3( 0 )
//	: mNumCtrlPoints( basis.mNumCtrlPoints ), mDegree( basis.mDegree ), mOpen( basis.mOpen ), mUniform( basis.mUniform )
{
    /*	int numKnots = initialize( mNumCtrlPoints, mDegree, mOpen );
     memcpy( mKnots, basis.mKnots, sizeof(float) * numKnots );*/
	*this = basis;
}

BSplineBasis& BSplineBasis::operator=( const BSplineBasis &basis )
{
    delete [] mKnots;
    deallocate2D( m_aafBD0 );
    deallocate2D( m_aafBD1 );
    deallocate2D( m_aafBD2 );
    deallocate2D( m_aafBD3 );
    
	mNumCtrlPoints = basis.mNumCtrlPoints;
	mDegree = basis.mDegree;
	mOpen = basis.mOpen;
	mUniform = basis.mUniform;
    
	if( mNumCtrlPoints > 0 ) {
		int numKnots = initialize( mNumCtrlPoints, mDegree, mOpen );
		memcpy( mKnots, basis.mKnots, sizeof(float) * numKnots );
	}
	else {
		mKnots = 0;
		m_aafBD0 = m_aafBD1 = m_aafBD2 = m_aafBD3 = 0;
	}
    
	return *this;
}

void BSplineBasis::create( int iNumCtrlPoints, int iDegree, bool bOpen )
{
    mUniform = true;
    
    int i, iNumKnots = initialize( iNumCtrlPoints, iDegree, bOpen );
    float fFactor = (1.0f)/( iNumCtrlPoints - mDegree );
    if ( mOpen ) {
        for ( i = 0; i <= mDegree; i++ ) {
            mKnots[i] = (float)0.0;
        }
        
        for (/**/; i < iNumCtrlPoints; i++ ) {
            mKnots[i] = ( i - mDegree ) * fFactor;
        }
        
        for(/**/; i < iNumKnots; i++) {
            mKnots[i] = (float)1.0;
        }
    }
    else {
        for ( i = 0; i < iNumKnots; i++ ) {
            mKnots[i] = ( i - mDegree ) * fFactor;
        }
    }
}

BSplineBasis::BSplineBasis( int aNumCtrlPoints, int iDegree, const float *afKnot )
{
    create( mNumCtrlPoints, iDegree, afKnot );
}

void BSplineBasis::create( int aNumCtrlPoints, int iDegree, const float *afKnot )
{
    mUniform = false;
    
	mNumCtrlPoints = aNumCtrlPoints;
    
    int i, iNumKnots = initialize( mNumCtrlPoints, iDegree, true );
    for( i = 0; i <= mDegree; i++ ) {
        mKnots[i] = (float)0.0;
    }
    
    for( int j = 0; i < mNumCtrlPoints; i++, j++ ) {
        mKnots[i] = afKnot[j];
    }
    
    for( /**/; i < iNumKnots; i++ ) {
        mKnots[i] = (float)1.0;
    }
}

BSplineBasis::~BSplineBasis()
{
    delete [] mKnots;
    deallocate2D( m_aafBD0 );
    deallocate2D( m_aafBD1 );
    deallocate2D( m_aafBD2 );
    deallocate2D( m_aafBD3 );
}

int BSplineBasis::getNumControlPoints() const
{
    return mNumCtrlPoints;
}

int BSplineBasis::getDegree() const
{
    return mDegree;
}


bool BSplineBasis::isOpen() const
{
    return mOpen;
}


bool BSplineBasis::isUniform() const
{
    return mUniform;
}


float BSplineBasis::getD0( int i ) const
{
    return m_aafBD0[mDegree][i];
}


float BSplineBasis::getD1( int i ) const
{
    return m_aafBD1[mDegree][i];
}


float BSplineBasis::getD2( int i ) const
{
    return m_aafBD2[mDegree][i];
}


float BSplineBasis::getD3( int i ) const
{
    return m_aafBD3[mDegree][i];
}


float** BSplineBasis::allocate() const
{
    int iRows = mDegree + 1;
    int iCols = mNumCtrlPoints + mDegree;
    float** aafArray;
    allocate2D<float>( iCols, iRows, aafArray );
    memset(aafArray[0],0,iRows*iCols*sizeof(float));
    return aafArray;
}

int BSplineBasis::initialize( int iNumCtrlPoints, int iDegree, bool bOpen )
{
    assert(iNumCtrlPoints >= 2);
    assert(1 <= iDegree && iDegree <= iNumCtrlPoints-1);
    
    mNumCtrlPoints = iNumCtrlPoints;
    mDegree = iDegree;
    mOpen = bOpen;
    
    int iNumKnots = mNumCtrlPoints+mDegree+1;
    mKnots = new float[iNumKnots];
    
    m_aafBD0 = allocate();
    m_aafBD1 = 0;
    m_aafBD2 = 0;
    m_aafBD3 = 0;
    
    return iNumKnots;
}

void BSplineBasis::setKnot( int i, float fKnot )
{
	if( ! mUniform ) {
		// access only allowed to elements d+1 <= j <= n
		int j = i + mDegree + 1;
		if( mDegree+1 <= j && j <= mNumCtrlPoints - 1 ) {
			mKnots[j] = fKnot;
		}
	}
}

float BSplineBasis::getKnot( int i ) const
{
    if( ! mUniform ) {
        // access only allowed to elements d+1 <= j <= n
        int j = i + mDegree + 1;
        if ( ( mDegree + 1 <= j ) && ( j <= mNumCtrlPoints - 1 ) ) {
            return mKnots[j];
        }
    }
    
    return std::numeric_limits<float>::max();
}


int BSplineBasis::getKey( float& rfTime ) const
{
	if( mOpen ) {
		// open splines clamp to [0,1]
		if( rfTime <= (float)0.0 ) {
			rfTime = (float)0.0;
			return mDegree;
		}
		else if ( rfTime >= (float)1.0 ) {
			rfTime = (float)1.0;
			return mNumCtrlPoints - 1;
		}
	}
	else {
		// periodic splines wrap to [0,1]
		if (rfTime < (float)0.0 || rfTime >= (float)1.0) {
			rfTime -= floorf( rfTime );
		}
	}
    
    
	int i;
	if( mUniform ) {
		i = mDegree + (int)( ( mNumCtrlPoints - mDegree ) * rfTime );
	}
	else {
		for( i = mDegree + 1; i <= mNumCtrlPoints; i++ ) {
			if( rfTime < mKnots[i] ) {
				break;
			}
		}
		i--;
	}
    
	return i;
}

void BSplineBasis::compute( float fTime, unsigned int uiOrder, int &riMinIndex, int &riMaxIndex ) const
{
    // only derivatives through third order currently supported
    assert(uiOrder <= 3);
    
    if (uiOrder >= 1) {
        if (!m_aafBD1) {
            m_aafBD1 = allocate();
        }
        
        if (uiOrder >= 2) {
            if( ! m_aafBD2 ) {
                m_aafBD2 = allocate();
            }
            
            if( uiOrder >= 3 ) {
                if ( ! m_aafBD3 ) {
                    m_aafBD3 = allocate();
                }
            }
        }
    }
    
    int i = getKey(fTime);
    m_aafBD0[0][i] = (float)1.0;
    
    if( uiOrder >= 1 ) {
        m_aafBD1[0][i] = (float)0.0;
        if ( uiOrder >= 2 ) {
            m_aafBD2[0][i] = (float)0.0;
            if ( uiOrder >= 3 ) {
                m_aafBD3[0][i] = (float)0.0;
            }
        }
    }
    
    float fN0 = fTime-mKnots[i], fN1 = mKnots[i+1] - fTime;
    float fInvD0, fInvD1;
    int j;
    for( j = 1; j <= mDegree; j++ ) {
		fInvD0 = ((float)1.0)/(mKnots[i+j]-mKnots[i]);
		fInvD1 = ((float)1.0)/(mKnots[i+1]-mKnots[i-j+1]);
        
		m_aafBD0[j][i] = fN0*m_aafBD0[j-1][i]*fInvD0;
		m_aafBD0[j][i-j] = fN1*m_aafBD0[j-1][i-j+1]*fInvD1;
        
		if( uiOrder >= 1 ) {
			m_aafBD1[j][i] = (fN0*m_aafBD1[j-1][i]+m_aafBD0[j-1][i])*fInvD0;
			m_aafBD1[j][i-j] = (fN1*m_aafBD1[j-1][i-j+1]-m_aafBD0[j-1][i-j+1])
            *fInvD1;
            
			if( uiOrder >= 2 ) {
				m_aafBD2[j][i] = (fN0*m_aafBD2[j-1][i] +
                                  ((float)2.0)*m_aafBD1[j-1][i])*fInvD0;
				m_aafBD2[j][i-j] = (fN1*m_aafBD2[j-1][i-j+1] -
                                    ((float)2.0)*m_aafBD1[j-1][i-j+1])*fInvD1;
                
				if ( uiOrder >= 3 ) {
					m_aafBD3[j][i] = (fN0*m_aafBD3[j-1][i] +
                                      ((float)3.0)*m_aafBD2[j-1][i])*fInvD0;
					m_aafBD3[j][i-j] = (fN1*m_aafBD3[j-1][i-j+1] -
                                        ((float)3.0)*m_aafBD2[j-1][i-j+1])*fInvD1;
				}
			}
		}
    }
    
    for( j = 2; j <= mDegree; j++ ) {
        for( int k = i-j+1; k < i; k++ ) {
            fN0 = fTime-mKnots[k];
            fN1 = mKnots[k+j+1]-fTime;
            fInvD0 = ((float)1.0)/(mKnots[k+j]-mKnots[k]);
            fInvD1 = ((float)1.0)/(mKnots[k+j+1]-mKnots[k+1]);
            
            m_aafBD0[j][k] = fN0*m_aafBD0[j-1][k]*fInvD0 + fN1*
            m_aafBD0[j-1][k+1]*fInvD1;
            
            if( uiOrder >= 1 ) {
				m_aafBD1[j][k] = (fN0*m_aafBD1[j-1][k]+m_aafBD0[j-1][k])*
                fInvD0 + (fN1*m_aafBD1[j-1][k+1]-m_aafBD0[j-1][k+1])*
                fInvD1;
                
				if( uiOrder >= 2 ) {
					m_aafBD2[j][k] = (fN0*m_aafBD2[j-1][k] +
                                      ((float)2.0)*m_aafBD1[j-1][k])*fInvD0 +
                    (fN1*m_aafBD2[j-1][k+1]- ((float)2.0)*
                     m_aafBD1[j-1][k+1])*fInvD1;
                    
                    if( uiOrder >= 3 ) {
						m_aafBD3[j][k] = (fN0*m_aafBD3[j-1][k] +
                                          ((float)3.0)*m_aafBD2[j-1][k])*fInvD0 +
                        (fN1*m_aafBD3[j-1][k+1] - ((float)3.0)*
                         m_aafBD2[j-1][k+1])*fInvD1;
                    }
                }
            }
        }
    }
    
    riMinIndex = i - mDegree;
    riMaxIndex = i;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Integration functions for arc length
template<typename T>
typename T::TYPE getSpeedWithData( float fTime, void* pvData)
{
    return ((BSpline<T>*)pvData)->getSpeed( fTime );
}

template<typename Real>
Real rombergIntegral( int iOrder, Real fA, Real fB, Real (*oF)(Real,void*), void* pvUserData )
{
    assert(iOrder > 0);
    Real** aafRom;
    //Allocate<Real>(iOrder,2,aafRom);
    allocate2D<Real>( iOrder, 2, aafRom );
    
    Real fH = fB - fA;
    
    aafRom[0][0] = ((Real)0.5)*fH*(oF(fA,pvUserData)+oF(fB,pvUserData));
    for (int i0=2, iP0=1; i0 <= iOrder; i0++, iP0 *= 2, fH *= (Real)0.5)
    {
        // approximations via the trapezoid rule
        Real fSum = (Real)0.0;
        int i1;
        for (i1 = 1; i1 <= iP0; i1++)
        {
            fSum += oF(fA + fH*(i1-((Real)0.5)),pvUserData);
        }
        
        // Richardson extrapolation
        aafRom[1][0] = ((Real)0.5)*(aafRom[0][0] + fH*fSum);
        for (int i2 = 1, iP2 = 4; i2 < i0; i2++, iP2 *= 4)
        {
            aafRom[1][i2] = (iP2*aafRom[1][i2-1] - aafRom[0][i2-1])/(iP2-1);
        }
        
        for (i1 = 0; i1 < i0; i1++)
        {
            aafRom[0][i1] = aafRom[1][i1];
        }
    }
    
    Real fResult = aafRom[0][iOrder-1];
    deallocate2D<Real>( aafRom );
    return fResult;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// BSpline
template<typename T>
BSpline<T>::BSpline( const std::vector<T> &points, int degree, bool loop, bool open )
: mLoop( loop )
{
	assert(points.size() >= 2);
	assert( ( 1 <= degree ) && ( degree <= (int)points.size() - 1 ) );
    
	mNumCtrlPoints = (int)points.size();
	mReplicate = ( mLoop ? (open ? 1 : degree) : 0);
	createControl( &points[0] );
	mBasis.create( mNumCtrlPoints + mReplicate, degree, open );
}

template<typename T>
BSpline<T>::BSpline( int numControlPoints, const T *controlPoints, int degree, bool loop, const float *knots )
: mNumCtrlPoints( numControlPoints ), mLoop( loop )
{
	assert( mNumCtrlPoints >= 2);
	assert( ( 1 <= degree ) && ( degree <= mNumCtrlPoints - 1 ) );
    
	mReplicate = (mLoop ? 1 : 0);
	createControl( controlPoints );
	mBasis.create( mNumCtrlPoints + mReplicate, degree, knots );
}

template<typename T>
BSpline<T>::BSpline( const BSpline &bspline )
: mCtrlPoints( 0 )
{
	*this = bspline;
}

template<typename T>
BSpline<T>& BSpline<T>::operator=( const BSpline &bspline )
{
	delete [] mCtrlPoints;
    
	mNumCtrlPoints = bspline.mNumCtrlPoints;
	mLoop = bspline.mLoop;
	mBasis = bspline.mBasis;
	mReplicate = bspline.mReplicate;
    
	if( mNumCtrlPoints > 0 )
		createControl( bspline.mCtrlPoints );
	else
		mCtrlPoints = 0;
    
	return *this;
}


template<typename T>
BSpline<T>::~BSpline()
{
	delete [] mCtrlPoints;
}

template<typename T>
void BSpline<T>::createControl( const T* akCtrlPoint )
{
	int iNewNumCtrlPoints = mNumCtrlPoints + mReplicate;
	mCtrlPoints = new T[iNewNumCtrlPoints];
	size_t uiSrcSize = mNumCtrlPoints*sizeof(T);
	memcpy( mCtrlPoints, akCtrlPoint, uiSrcSize );
	for( int i = 0; i < mReplicate; i++ ) {
		mCtrlPoints[mNumCtrlPoints+i] = akCtrlPoint[i];
	}
}

template<typename T>
void BSpline<T>::setControlPoint( int i, const T& rkCtrl )
{
	if( ( 0 <= i ) && ( i < mNumCtrlPoints ) ) {
		// set the control point
		mCtrlPoints[i] = rkCtrl;
        
		// set the replicated control point
		if( i < mReplicate ) {
			mCtrlPoints[mNumCtrlPoints+i] = rkCtrl;
		}
	}
}

template<typename T>
T BSpline<T>::getControlPoint( int i ) const
{
    if( ( 0 <= i ) && ( i < mNumCtrlPoints) ) {
		return mCtrlPoints[i];
    }
    
    return mCtrlPoints[0]; // zach
}

template<typename T>
void BSpline<T>::setKnot( int i, float fKnot )
{
    mBasis.setKnot( i, fKnot );
}

template<typename T>
float BSpline<T>::getKnot( int i ) const
{
    return mBasis.getKnot( i );
}

template<typename T>
void BSpline<T>::get( float t, T *position, T *firstDerivative, T *secondDerivative, T *thirdDerivative ) const
{
	int i, iMin, iMax;
	if( thirdDerivative ) {
		mBasis.compute( t, 3, iMin, iMax );
	}
	else if( secondDerivative ) {
		mBasis.compute( t, 2, iMin, iMax );
	}
	else if( firstDerivative ) {
		mBasis.compute( t, 1, iMin, iMax );
	}
	else {
        mBasis.compute( t, 0, iMin, iMax );
	}
    
	if( position ) {
		*position = T::zero();
		for( i = iMin; i <= iMax; i++ ) {
			float weight = mBasis.getD0( i );
			*position += mCtrlPoints[i] * weight;
		}
	}
    
	if( firstDerivative ) {
		*firstDerivative = T::zero();
		for( i = iMin; i <= iMax; i++ ) {
			*firstDerivative += mCtrlPoints[i]*mBasis.getD1( i );
		}
	}
    
	if( secondDerivative ) {
		*secondDerivative = T::zero();
		for( i = iMin; i <= iMax; i++ ) {
			*secondDerivative += mCtrlPoints[i]*mBasis.getD2( i );
		}
	}
    
	if( thirdDerivative ) {
		*thirdDerivative = T::zero();
		for (i = iMin; i <= iMax; i++) {
			*thirdDerivative += mCtrlPoints[i]*mBasis.getD3( i );
		}
	}
}

template<typename T>
float BSpline<T>::getTime( float length ) const
{
	const size_t MAX_ITERATIONS = 32;
	const float TOLERANCE = 1.0e-03f;
	// ensure that we remain within valid parameter space
	float totalLength = getLength( 0, 1 );
	if( length >= totalLength )
		return 1;
	if( length <= 0 )
		return 0;
    
	// initialize bisection endpoints
	float a = 0, b = 1;
	float p = length / totalLength;    // make first guess
    
	// iterate and look for zeros
	for ( size_t i = 0; i < MAX_ITERATIONS; ++i ) {
		// compute function value and test against zero
		float func = getLength( 0, p ) - length;
		if(abs( func ) < TOLERANCE ) {
			return p;
		}
        
        // update bisection endpoints
		if( func < 0 ) {
			a = p;
		}
		else {
			b = p;
		}
        
		// get speed along curve
		float speed = 1.0;// getSpeed( p );  // zach ??
        
		// if result will lie outside [a,b]
		if( ((p-a)*speed - func)*((p-b)*speed - func) > -TOLERANCE ) {
			// do bisection
			p = 0.5f*(a+b);
		}
		else {
			// otherwise Newton-Raphson
			p -= func/speed;
		}
	}
    
	// We failed to converge, but hopefully 'p' is close enough anyway
	return p;
}

template<typename T>
float BSpline<T>::getLength( float fT0, float fT1 ) const
{
	if( fT0 >= fT1 )
		return (float)0.0;
    
    return 0; //rombergIntegral<typename T::TYPE>( 10, fT0, fT1, getSpeedWithData<T>, (void*)this );
}

template<typename T>
BSplineBasis& BSpline<T>::getBasis()
{
	return mBasis;
}

template<typename T>
T BSpline<T>::getPosition( float t ) const
{
	T kPos;
	get( t, &kPos, 0, 0, 0 );
	return kPos;
}

template<typename T>
T BSpline<T>::getDerivative( float t ) const
{
	T kDer1;
	get( t, 0, &kDer1, 0, 0 );
	return kDer1;
}

template<typename T>
T BSpline<T>::getSecondDerivative( float t ) const
{
	T kDer2;
	get( t, 0, 0, &kDer2, 0 );
	return kDer2;
}

template<typename T>
T BSpline<T>::getThirdDerivative( float t ) const
{
	T kDer3;
	get( t, 0, 0, 0, &kDer3 );
	return kDer3;
}

// explicit template instantiations
template class BSpline<ofVec2f>;
template class BSpline<ofVec3f>;
template class BSpline<ofVec4f>;


    
    template<typename T>
    class BSplineFit
    {
    public:
        // Construction and destruction.  The preconditions for calling the
        // constructor are
        //   1 <= iDegree && iDegree < iControlQuantity <= iSampleQuantity
        // The samples point are contiguous blocks of iDimension real value
        // stored in afSampleData.
        BSplineFit( int iDimension, int iSampleQuantity, const T* afSampleData, int iDegree, int iControlQuantity );
        ~BSplineFit();
        
        // Access to input sample information.
        int getDimension() const;
        int getSampleQuantity() const;
        const T* getSampleData() const;
        
        // Access to output control point and curve information.
        int getDegree() const;
        int getControlQuantity() const;
        const T* getControlData() const;
        const BSplineFitBasis<T>& getBasis() const;
        
        // Evaluation of the B-spline curve.  It is defined for 0 <= t <= 1.  If
        // a t-value is outside [0,1], an open spline clamps it to [0,1].  The
        // caller must ensure that afPosition[] has (at least) 'dimension'
        // elements.
        void getPosition( T fT, T *afPosition ) const;
        
    private:
        // The matric inversion calculations are performed with double-precision,
        // even when the type T is 'float'.
        bool choleskyFactor( BandedMatrixd &rkMatrix ) const;
        bool solveLower( BandedMatrixd &rkMatrix, double* adControlData ) const;
        bool solveUpper( BandedMatrixd &rkMatrix, double* adControlData ) const;
        
        // Input sample information.
        int m_iDimension;
        int m_iSampleQuantity;
        const T* m_afSampleData;
        
        // The fitted B-spline curve, open and with uniform knots.
        int m_iDegree;
        int m_iControlQuantity;
        T* m_afControlData;
        BSplineFitBasis<T> m_kBasis;
    };
    
    typedef BSplineFit<float> BSplineFitf;
    typedef BSplineFit<double> BSplineFitd;
    typedef BSplineFitBasis<float> BSplineFitBasisf;
    typedef BSplineFitBasis<double> BSplineFitBasisd;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    template<typename T>
    BSplineFitBasis<T>::BSplineFitBasis(int iQuantity, int iDegree)
    {
        assert(1 <= iDegree && iDegree < iQuantity);
        m_iQuantity = iQuantity;
        m_iDegree = iDegree;
        
        m_afValue = new T[iDegree+1];
        m_afKnot = new T[2*iDegree];
    }
    
    template<typename T>
    BSplineFitBasis<T>::~BSplineFitBasis()
    {
        delete [] m_afValue;
        delete [] m_afKnot;
    }
    
    template<typename T>
    int BSplineFitBasis<T>::getQuantity() const
    {
        return m_iQuantity;
    }
    
    template<typename T>
    int BSplineFitBasis<T>::getDegree() const
    {
        return m_iDegree;
    }
    
    template<typename T>
    void BSplineFitBasis<T>::compute( T fTime, int &iMin, int &iMax ) const
    {
        assert((T)0.0 <= fTime && fTime <= (T)1.0);
        
        // Use scaled time and scaled knots so that 1/(Q-D) does not need to
        // be explicitly stored by the class object.  Determine the extreme
        // indices affected by local control.
        T fQmD = (T)(m_iQuantity - m_iDegree);
        T fT;
        if (fTime <= (T)0.0)
        {
            fT = (T)0.0;
            iMin = 0;
            iMax = m_iDegree;
        }
        else if (fTime >= (T)1.0)
        {
            fT = fQmD;
            iMax = m_iQuantity - 1;
            iMin = iMax - m_iDegree;
        }
        else
        {
            fT = fQmD*fTime;
            iMin = (int)fT;
            iMax = iMin + m_iDegree;
        }
        
        // Precompute the knots.
        for (int i0 = 0, i1 = iMax+1-m_iDegree; i0 < 2*m_iDegree; i0++, i1++)
        {
            if (i1 <= m_iDegree)
            {
                m_afKnot[i0] = (T)0.0f;
            }
            else if (i1 >= m_iQuantity)
            {
                m_afKnot[i0] = fQmD;
            }
            else
            {
                m_afKnot[i0] = (T)(i1 - m_iDegree);
            }
        }
        
        // Initialize the basis function evaluation table.  The first degree-1
        // entries are zero, but they do not have to be set explicitly.
        m_afValue[m_iDegree] = (T)1.0;
        
        // Update the basis function evaluation table, each iteration overwriting
        // the results from the previous iteration.
        for (int iRow = m_iDegree-1; iRow >= 0; iRow--)
        {
            int iK0 = m_iDegree, iK1 = iRow;
            T fKnot0 = m_afKnot[iK0], fKnot1 = m_afKnot[iK1];
            T fInvDenom = ((T)1.0)/(fKnot0 - fKnot1);
            T fC1 = (fKnot0 - fT)*fInvDenom, fC0;
            m_afValue[iRow] = fC1*m_afValue[iRow+1];
            
            for (int iCol = iRow+1; iCol < m_iDegree; iCol++)
            {
                fC0 = (fT - fKnot1)*fInvDenom;
                m_afValue[iCol] *= fC0;
                
                fKnot0 = m_afKnot[++iK0];
                fKnot1 = m_afKnot[++iK1];
                fInvDenom = ((T)1.0)/(fKnot0 - fKnot1);
                fC1 = (fKnot0 - fT)*fInvDenom;
                m_afValue[iCol] += fC1*m_afValue[iCol+1];
            }
            
            fC0 = (fT - fKnot1)*fInvDenom;
            m_afValue[m_iDegree] *= fC0;
        }
    }
    
    template<typename T>
    T BSplineFitBasis<T>::getValue( int i ) const
    {
        assert(0 <= i && i <= m_iDegree);
        return m_afValue[i];
    }
    
    //----------------------------------------------------------------------------
    
    template<typename T>
    BSplineFit<T>::BSplineFit( int iDimension, int iSampleQuantity, const T* afSampleData, int iDegree, int iControlQuantity )
    : m_kBasis( iControlQuantity, iDegree )
    {
        if( iControlQuantity <= iDegree + 1 ) iControlQuantity = iDegree + 2;
        if( iControlQuantity > iSampleQuantity ) iControlQuantity = iSampleQuantity;
        iDegree = ofClamp( iDegree, 1, iControlQuantity - 1 );
        
        assert(iDimension >= 1);
        assert(1 <= iDegree && iDegree < iControlQuantity);
        assert(iControlQuantity <= iSampleQuantity);
        
        m_iDimension = iDimension;
        m_iSampleQuantity = iSampleQuantity;
        m_afSampleData = afSampleData;
        m_iDegree = iDegree;
        m_iControlQuantity = iControlQuantity;
        m_afControlData = new T[m_iDimension*iControlQuantity];
        
        // Fit the data points with a B-spline curve using a least-squares error
        // metric.  The problem is of the form A^T*A*X = A^T*B.
        BSplineFitBasisd kDBasis(m_iControlQuantity,m_iDegree);
        double dTMultiplier = 1.0/(double)(m_iSampleQuantity - 1);
        double dT;
        int i0, i1, i2, iMin, iMax, j;
        
        // Construct the matrix A (depends only on the output basis function).
        BandedMatrixd* pkAMat = new BandedMatrixd( m_iControlQuantity, m_iDegree+1, m_iDegree + 1 );
        
        for (i0 = 0; i0 < m_iControlQuantity; i0++) {
            for ( i1 = 0; i1 < i0; i1++ ) {
                (*pkAMat)(i0,i1) = (*pkAMat)(i1,i0);
            }
            
            int i1Max = i0 + m_iDegree;
            if (i1Max >= m_iControlQuantity)
            {
                i1Max = m_iControlQuantity - 1;
            }
            
            for (i1 = i0; i1 <= i1Max; i1++)
            {
                double dValue = 0.0;
                for (i2 = 0; i2 < m_iSampleQuantity; i2++)
                {
                    dT = dTMultiplier*(double)i2;
                    kDBasis.compute(dT,iMin,iMax);
                    if (iMin <= i0 && i0 <= iMax && iMin <= i1 && i1 <= iMax)
                    {
                        double dB0 = kDBasis.getValue(i0 - iMin);
                        double dB1 = kDBasis.getValue(i1 - iMin);
                        dValue += dB0*dB1;
                    }
                }
                (*pkAMat)(i0,i1) = dValue;
            }
        }
        
        // Construct the matrix B.
        double** aadBMat;
        aadBMat = new double*[m_iControlQuantity];
        aadBMat[0] = new double[m_iControlQuantity*m_iSampleQuantity];
        for( int iRow = 1; iRow < m_iControlQuantity; iRow++ ) {
            aadBMat[iRow] = &aadBMat[0][m_iSampleQuantity*iRow];
        }
        
        memset( aadBMat[0],0,m_iControlQuantity*m_iSampleQuantity*sizeof(double) );
        for (i0 = 0; i0 < m_iControlQuantity; i0++)
        {
            for (i1 = 0; i1 < m_iSampleQuantity; i1++)
            {
                dT = dTMultiplier*(double)i1;
                kDBasis.compute(dT,iMin,iMax);
                if (iMin <= i0 && i0 <= iMax)
                {
                    aadBMat[i0][i1] = kDBasis.getValue(i0 - iMin);
                }
            }
        }
        
        // Construct the control points for the least-squares curve.
        double* adControlData = new double[m_iDimension*m_iControlQuantity];
        memset( adControlData,0,m_iDimension*m_iControlQuantity*sizeof(double) );
        double* pdBaseTarget = adControlData;
        
        
        for (i0 = 0; i0 < m_iControlQuantity; i0++)
        {
            const T* pfSource = m_afSampleData;
            double* adTarget = pdBaseTarget;
            for (i1 = 0; i1 < m_iSampleQuantity; i1++)
            {
                double dBValue = aadBMat[i0][i1];
                for (j = 0; j < m_iDimension; j++)
                {
                    adTarget[j] += dBValue*(double)(*pfSource++);
                }
            }
            pdBaseTarget += m_iDimension;
        }
        
       
        // Solve A^T*A*ControlData = A^T*B*SampleData.
        bool bSolved = choleskyFactor(*pkAMat);
        assert(bSolved);
        bSolved = solveLower(*pkAMat,adControlData);
        assert(bSolved);
        bSolved = solveUpper(*pkAMat,adControlData);
        assert(bSolved);
        
        // Set the B-spline control points.
        T* pfTarget = m_afControlData;
        const double* pdSource = adControlData;
        for (i0 = 0; i0 < m_iDimension*m_iControlQuantity; i0++)
        {
            *pfTarget++ = (T)(*pdSource++);
        }
        
        // Set the first and last output control points to match the first and
        // last input samples.  This supports the application of fitting keyframe
        // data with B-spline curves.  The user expects that the curve passes
        // through the first and last positions in order to support matching two
        // consecutive keyframe sequences.
        T* pfCEnd0 = m_afControlData;
        const T* pfSEnd0 = m_afSampleData;
        T* pfCEnd1 = &m_afControlData[m_iDimension*(m_iControlQuantity-1)];
        const T* pfSEnd1 = &m_afSampleData[m_iDimension*(m_iSampleQuantity-1)];
        for (j = 0; j < m_iDimension; j++)
        {
            *pfCEnd0++ = *pfSEnd0++;
            *pfCEnd1++ = *pfSEnd1++;
        }
        
        delete [] adControlData;
        if( aadBMat ) {
            delete [] aadBMat[0];
            delete [] aadBMat;
            aadBMat = 0;
        }
        delete pkAMat;
    }
    
    template<typename T>
    BSplineFit<T>::~BSplineFit ()
    {
        delete [] m_afControlData;
    }
    
    template<typename T>
    int BSplineFit<T>::getDimension() const
    {
        return m_iDimension;
    }
    
    template<typename T>
    int BSplineFit<T>::getSampleQuantity() const
    {
        return m_iSampleQuantity;
    }
    
    template<typename T>
    const T* BSplineFit<T>::getSampleData() const
    {
        return m_afSampleData;
    }
    
    template<typename T>
    int BSplineFit<T>::getDegree() const
    {
        return m_iDegree;
    }
    
    template<typename T>
    int BSplineFit<T>::getControlQuantity() const
    {
        return m_iControlQuantity;
    }
    
    template<typename T>
    const T* BSplineFit<T>::getControlData() const
    {
        return m_afControlData;
    }
    
    template<typename T>
    const BSplineFitBasis<T>& BSplineFit<T>::getBasis() const
    {
        return m_kBasis;
    }
    
    template<typename T>
    void BSplineFit<T>::getPosition( T fT, T* afPosition ) const
    {
        assert(afPosition);
        
        int iMin, iMax;
        m_kBasis.compute(fT,iMin,iMax);
        
        T* pfSource = &m_afControlData[m_iDimension*iMin];
        T fBasisValue = m_kBasis.getValue(0);
        int j;
        for (j = 0; j < m_iDimension; j++)
        {
            afPosition[j] = fBasisValue*(*pfSource++);
        }
        
        for (int i = iMin+1, iIndex = 1; i <= iMax; i++, iIndex++)
        {
            fBasisValue = m_kBasis.getValue(iIndex);
            for (j = 0; j < m_iDimension; j++)
            {
                afPosition[j] += fBasisValue*(*pfSource++);
            }
        }
    }
    
    template<typename T>
    bool BSplineFit<T>::choleskyFactor( BandedMatrixd &rkMatrix ) const
    {
        const int iSize = rkMatrix.getSize(), iSizeM1 = iSize - 1;
        const int iBands = rkMatrix.getLBands();  // == GetUBands()
        
        int k, kMax;
        for (int i = 0; i < iSize; i++)
        {
            int jMin = i - iBands;
            if (jMin < 0)
            {
                jMin = 0;
            }
            
            int j;
            for (j = jMin; j < i; j++)
            {
                kMax = j + iBands;
                if (kMax > iSizeM1)
                {
                    kMax = iSizeM1;
                }
                
                for (k = i; k <= kMax; k++)
                {
                    rkMatrix(k,i) -= rkMatrix(i,j)*rkMatrix(k,j);
                }
            }
            
            kMax = j + iBands;
            if (kMax > iSizeM1)
            {
                kMax = iSizeM1;
            }
            
            for (k = 0; k < i; k++)
            {
                rkMatrix(k,i) = rkMatrix(i,k);
            }
            
            double dDiagonal = rkMatrix(i,i);
            if (dDiagonal <= 0.0)
            {
                return false;
            }
            double dInvSqrt = 1.0 / sqrt( dDiagonal );
            for (k = i; k <= kMax; k++)
            {
                rkMatrix(k,i) *= dInvSqrt;
            }
        }
        
        return true;
    }
    
    template<typename T>
    bool BSplineFit<T>::solveLower( BandedMatrixd& rkMatrix, double* adControlData ) const
    {
        const int iSize = rkMatrix.getSize();
        double* pdBaseTarget = adControlData;
        for (int iRow = 0; iRow < iSize; iRow++)
        {
            if( abs(rkMatrix(iRow,iRow)) < EPSILON_VALUE )
            {
                return false;
            }
            
            const double* pdBaseSource = adControlData;
            double* adTarget = pdBaseTarget;
            int j;
            for (int iCol = 0; iCol < iRow; iCol++)
            {
                const double* pdSource = pdBaseSource;
                double dMatValue = rkMatrix(iRow,iCol);
                for (j = 0; j < m_iDimension; j++)
                {
                    adTarget[j] -= dMatValue*(*pdSource++);
                }
                pdBaseSource += m_iDimension;
            }
            
            double dInverse = 1.0/rkMatrix(iRow,iRow);
            for (j = 0; j < m_iDimension; j++)
            {
                adTarget[j] *= dInverse;
            }
            pdBaseTarget += m_iDimension;
        }
        
        return true;
    }
    
    template<typename T>
    bool BSplineFit<T>::solveUpper( BandedMatrixd& rkMatrix, double *adControlData ) const
    {
        const int iSize = rkMatrix.getSize();
        double* pdBaseTarget = &adControlData[m_iDimension*(iSize-1)];
        for (int iRow = iSize - 1; iRow >= 0; iRow--)
        {
            if( abs(rkMatrix(iRow,iRow)) < EPSILON_VALUE ) {
                return false;
            }
            
            const double* pdBaseSource = &adControlData[m_iDimension*(iRow+1)];
            double* adTarget = pdBaseTarget;
            int j;
            for (int iCol = iRow+1; iCol < iSize; iCol++)
            {
                const double* pdSource = pdBaseSource;
                double dMatValue = rkMatrix(iRow,iCol);
                for (j = 0; j < m_iDimension; j++)
                {
                    adTarget[j] -= dMatValue*(*pdSource++);
                }
                pdBaseSource += m_iDimension;
            }
            
            double dInverse = 1.0/rkMatrix(iRow,iRow);
            for (j = 0; j < m_iDimension; j++)
            {
                adTarget[j] *= dInverse;
            }
            pdBaseTarget -= m_iDimension;
        }
        
        return true;
    }
    

    BSpline3f fitBSpline( std::vector<ofVec3f> &samples, int degree, int outputSamples )
    {
        
       // float * f = &(samples[0].x);
        BSplineFit<float> fit( 3, (int)samples.size(), samples[0].getPtr(), degree, outputSamples );
        
        vector<ofVec3f> points;
        for( int c = 0; c < fit.getControlQuantity(); ++c ) {
            
            ofVec3f vec;
            vec.set(fit.getControlData()[c * 3 + 0], fit.getControlData()[c * 3 + 1], fit.getControlData()[c * 3 + 2]);
            
            points.push_back( vec );
        }
        return BSpline3f( points, fit.getDegree(), false, true );
    }
    
    template class BSplineFit<float>;
    template class BSplineFit<double>;
    template class BSplineFitBasis<float>;
    template class BSplineFitBasis<double>;
    
//    template BSpline<ofVec2f> fitBSpline( const std::vector<ofVec2f> &samples, int degree, int outputSamples );
//    template BSpline<ofVec3f> fitBSpline( const std::vector<ofVec3f> &samples, int degree, int outputSamples );
//    template BSpline<Vec4f> fitBSpline( const std::vector<Vec4f> &samples, int degree, int outputSamples );
//    




//--------------------------------------------------------------
void testApp::setup(){

}

//--------------------------------------------------------------
void testApp::update(){

    
    
}

//--------------------------------------------------------------
void testApp::draw(){

    ofSetColor(255);
    //temp.draw();
    
    ofSetColor(255,0,0);
    if (temp.getVertices().size() > 15){
        BSpline3f spline = fitBSpline( temp.getVertices(), 3, temp.getVertices().size() *0.8);
        
        cout << spline.getNumControlPoints() << endl;
        ofMesh mesh;
        int count =  ofMap(mouseX, 0, ofGetWidth(), 3, 10000, true);
        mesh.setMode(OF_PRIMITIVE_LINE_STRIP);
        for (int i = 0; i < count; i++){
            mesh.addVertex(spline.getPosition( i / (float)count));
        }
        mesh.draw();
    }
    
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){

}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

    
    temp.addVertex( ofPoint(x,y));
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    temp.clear();
    temp.addVertex( ofPoint(x,y));
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}
