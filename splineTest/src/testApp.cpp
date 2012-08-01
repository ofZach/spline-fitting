#include "testApp.h"
#include "BSpline.h"




const double BoundaryConditions[3][4] =
{
    //  0       1       M-1     M
    {
        -4,
        -1,
        -1,
        -4 },
    {
        0,
        1,
        1,
        0 },
    {
        2,
        -1,
        -1,
        2 } };


inline double Beta(int m, int M)
{
    if (m > 1 && m < M-1)
        return 0.0;
    if (m >= M-1)
        m -= M-3;
    return BoundaryConditions[0][m];
}

double Basis(int m, int M, float x, float xmin, float DX){
    double y = 0;
    double xm = xmin + (m * DX);
    double z = abs((double)(x - xm) / (double)DX);
    if (z < 2.0) {
        z = 2 - z;
        y = 0.25 * (z*z*z);
        z -= 1.0;
        if (z > 0)
            y -= (z*z*z);
    }
    
    // Boundary conditions, if any, are an additional addend.
    if (m == 0 || m == 1)
        y += Beta(m, M) * Basis(-1, M, x, xmin, DX);
    else if (m == M-1 || m == M)
        y += Beta(m, M) * Basis(M+1, M, x, xmin, DX);
    
    return y;
}

float evaluate(float x, float * coeficients, float xmin, float DX, int M, float mean) {
    float y = 0;
    if (true) {
        int n = (int)((x - xmin)/DX);
        
		// M ? 
		// s-A[] ?
		
        for (int i = max(0, n-1); i <= min(M, n+2); ++i) {
            y += coeficients[i] * Basis(i, M, x, xmin, DX);
        }
        y += mean;
    }
    return y;
}


//--------------------------------------------------------------
//--------------------------------------------------------------



typedef double datum;
typedef BSpline<datum> SplineT;

typedef struct {
    
    int M;
    float xmin;
    float xmax;
    float DX;
    float mean;
    vector < float > coefficients;
    
} splineInfo;


vector < ofPoint > pts;
vector < ofPoint> ptsResult;
vector<datum> xVals;
vector<datum> yVals;


//--------------------------------------------------------------
//--------------------------------------------------------------
void testApp::setup(){

    
    
    for (int i = 0; i < 100; i++){
        
        float x = i;
        float y = ofNoise(1000 + i / 10.0) * 50;
        
        pts.push_back(ofPoint(x,y));
        
        xVals.push_back(x);
        yVals.push_back(y);
        
    }
    
    /*
     
     
     SplineT::Debug(0);
     SplineT spline(&xVals[0],
     xVals.size(),
     &yVals[0],
     0,
     0,
     12);
     
     printf("coefficients \n");
     
     splineInfo info;
     
     for (int i = 0; i <= spline.M; i++){
     printf("%i = %f \n", i, spline.coefficient(i));
     info.coefficients.push_back(spline.coefficient(i));
     info.M = spline.M;
     info.DX = spline.DX;
     info.xmin = spline.xmin;
     info.xmax = spline.xmax;
     info.mean = spline.mean;
     }
     
     
     for (int i = 0; i < 100; i++){
     
     float y = evaluate(i, &info.coefficients[0], info.xmin, info.DX, info.M, info.mean);
     //evaluate(<#float x#>, <#float *coeficients#>, <#float xmin#>, <#float DX#>, <#int M#>, <#float mean#>)
     ptsResult.push_back(ofPoint(i,y));
     
     }
     
     
     
     */
    
    

    
    
}

//--------------------------------------------------------------
void testApp::update(){

}

//--------------------------------------------------------------
void testApp::draw(){

    
    ofPushMatrix();
    
    ofTranslate(50,50);
    
    ofNoFill();
    ofBeginShape();
    for (int i = 0; i < pts.size(); i++){
        ofVertex(pts[i].x, pts[i].y);
    }
    ofEndShape();
    
    ofPopMatrix();
    
    
    
    
    
    
    ofTranslate(250,50);
    
    for (int j = 5; j < 50; j+=3){

        SplineT::Debug(0);
        SplineT spline(&xVals[0],
                       xVals.size(),
                       &yVals[0],
                       0,
                       0,
                       j);
        
        //printf("coefficients \n");
        
        splineInfo info;
        
        for (int i = 0; i <= spline.M; i++){
            //printf("%i = %f \n", i, spline.coefficient(i));
            info.coefficients.push_back(spline.coefficient(i));
            info.M = spline.M;
            info.DX = spline.DX;
            info.xmin = spline.xmin;
            info.xmax = spline.xmax;
            info.mean = spline.mean;
        }
        
        
        ofBeginShape();
        for (int i = 0; i < 100; i++){
            
            float y = evaluate(i, &info.coefficients[0], info.xmin, info.DX, info.M, 0); //info.mean);   // not really sure about the mean ?
            ptsResult.push_back(ofPoint(i,y));
            ofVertex(i,y);
        }
        ofEndShape();
        
        for (int i = 0; i < info.coefficients.size(); i++){
            
            float x = (i * 100.0/(info.coefficients.size()-1));
            float y = info.coefficients[i] * mouseX;
            //ofRect(x,y,2,2);
        }
        
        ofDrawBitmapStringHighlight("nPts = " + ofToString(j), ofPoint(150,0));
        
        ofTranslate(0,50);
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

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

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