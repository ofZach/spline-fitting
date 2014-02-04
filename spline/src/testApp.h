#pragma once

#include "ofMain.h"

template<typename T>
class BandedMatrix {
public:
	BandedMatrix( int iSize, int iLBands, int iUBands );
	BandedMatrix( const BandedMatrix& rkM );
	~BandedMatrix();
    
	BandedMatrix& operator=( const BandedMatrix &rkM );
    
	int getSize() const;
	int getLBands() const;
	int getUBands() const;
    
	T* getDBand();
	const T* getDBand() const;
    
	int getLBandMax( int i ) const;  // LBand(i):  0 <= index < LBandMax
	T* getLBand( int i );
	const T* getLBand( int i ) const;
    
	int getUBandMax( int i ) const;  // UBand(i):  0 <= index < UBandMax
	T* getUBand( int i );
	const T* getUBand( int i ) const;
    
	T& operator()( int iRow, int iCol );
	T operator()( int iRow, int iCol ) const;
    
	void setZero();
	void setIdentity();
    
private:
	void allocate();
	void deallocate();
    
	int m_iSize, m_iLBands, m_iUBands;
	T* m_afDBand;
	T** m_aafLBand;
	T** m_aafUBand;
};

typedef BandedMatrix<float> BandedMatrixf;
typedef BandedMatrix<double> BandedMatrixd;


class testApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
    
    
    
        ofPolyline temp;
    
		
};
