

//    points[0] = Point3Df(20,20,0) ;
//    points[1] = Point3Df(20,80,0) ;
//    points[2] = Point3Df(20,120,0) ;
//    points[3] = Point3Df(20,160,0) ;
//    points[4] = Point3Df(80,200,0) ;
//    points[5] = Point3Df(120,200,0) ;
//    points[6] = Point3Df(160,160,0) ;
//    points[7] = Point3Df(160,120,0) ;
//    points[8] = Point3Df(120,80,0) ;
//    points[9] = Point3Df(80,80,0) ;
//    
#include "nurbs.h"
#include "testApp.h"
#include "nurbs.cpp"


PlNurbsCurvef curveA,curveB,curveC,curveD;
Vector_Point3Df points(100) ;
using namespace PLib ; 
vector < ofPoint > pts;

//--------------------------------------------------------------
//--------------------------------------------------------------
void testApp::setup(){
    
    
    
    
    for (int i = 0; i < 100; i++){
        
        float x = i;
        float y = ofNoise(2000 + i / 20.0) * 50;
        
        points[i] = Point3Df(x,y,0);
        
        pts.push_back(ofPoint(x,y,0));
        
        
    }
    
    
    //curveA.globalInterp(points,3) ;
    //curveB.leastSquares(points,3,5) ;
    //curveC.leastSquares(points,3,8) ;
    
    
    
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
        
        curveB.leastSquares(points,3,j);
        
        ofNoFill();
        ofBeginShape();
        for (int i = 0; i < 100; i++){
            float pct = (float)i/100.0f;
            float x = curveB.pointAt(pct).data[0];
            float y = curveB.pointAt(pct).data[1];
            cout << x << " " << y << endl;
            ofVertex(x, y);
        }
        ofEndShape();
        ofFill();
        
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