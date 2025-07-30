#include "ofApp.h"



extern "C" {
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include "DataExtraction.h"
#include "SolarSystem.h"
}





void ofApp::coldStartInitialConditionzzz(Body* &simBodies, size_t numBodies)
{
	Vec3D position, velocity;
	double mass = 1.0;
	for (size_t i = 0; i < numBodies; i++)
	{
		//position = Vec3D((0 + rand() % (1500)), (0 + rand() % (1500)), (0 + rand() % (1500)));
		//velocity = Vec3D((0 + rand() % (1)), (0 + rand() % (1)), (0 + rand() % (1)));
		
		position = Vec3D(ofRandom(-1500, 1500), ofRandom(-1500, 1500), ofRandom(-1500, 1500));
		velocity = Vec3D(ofRandom(-1, 1), ofRandom(-1, 1), ofRandom(-1, 1));
		
		
		simBodies[i] = {position, velocity, mass};
	}
	
}


//--------------------------------------------------------------
void ofApp::setup()
{
	ofSetFrameRate(120);
	ofSetBackgroundColor(11);
	ofSetCircleResolution(32);
	ofEnableAntiAliasing();
	ofSetEscapeQuitsApp(false);
	ofSetVerticalSync(true);
	ofEnableSmoothing();
	ofEnableDepthTest();
	//ofDisableSmoothing();
	//ofDisableAntiAliasing();
	
	///*
	 int count = 0;
	 EphemerisRecord* marsData = parseHorizonsFilez("/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_mars.txt", &count);
	 if(!marsData)
	 {
	 printf("No records parsed.\n");
	 }
	 
	 printf("Parsed %d records.\n", count);
	 if(count > 0) {
	 printf("First record:\n");
	 printf("  JD=%.6f\n", marsData[0].julianDate);
	 printf("  Position (km): %g, %g, %g\n",
	 marsData[0].position[0],
	 marsData[0].position[1],
	 marsData[0].position[2]);
	 printf("  Velocity (km/s): %g, %g, %g\n",
	 marsData[0].velocity[0],
	 marsData[0].velocity[1],
	 marsData[0].velocity[2]);
	 }
	 
	 // Convert to meters if you wish:
	 // for(int i=0; i<count; i++){
	 //     marsData[i].position[0] *= 1000.0;
	 //     marsData[i].position[1] *= 1000.0;
	 //     ...
	 // }
	 
	 free(marsData);
	 //*/
	
	numBodies = 4444;
	bodies = new Body[numBodies];
	
	
	
	//Vec3D galaxyAnchorPoint = {0.0 + ofGetScreenWidth()*0.5, 0.0 + ofGetScreenHeight()*0.5, 0};
	//plummerModelInitialConditions(bodies, numBodies, 100000, 100000, 1, 1500, Vec3D(ofRandom(-5000, 5000), ofRandom(-5000, 5000), ofRandom(-5000, 5000))); // galaxyAnchorPoint);
	//plummerModelInitialConditions(bodies, numBodies, numBodies*4, numBodies*4, 1, 150, Vec3D(ofRandom(-150, 150), ofRandom(-150, 150), ofRandom(-150, 150))); // galaxyAnchorPoint);
	
	
	//coldStartInitialConditionzzz(bodies, numBodies);
	//coldStartInitialConditions(bodies, numBodies);
	//ColdStart(bodies, numBodies);
	SetupSolarSystem(bodies, numBodies);
	
	
	//ColdStart(bodies, numBodies);
	
	//
	/*
	 for (size_t i = 0; i < numBodies; i++)
	 {
	 bodies[i] = Body(Vec3D(ofRandom(1500), ofRandom(1500), ofRandom(1500)), Vec3D(ofRandom(1), ofRandom(1), ofRandom(1)), 1);
	 }
	 //*/
	
	
	//
	/*
	 for (size_t i = 0; i < numBodies; i++)
	 {
	 bodies[i] = Body(Vec3D((arc4random()%(3000))-1500, (arc4random()%(3000))-1500, (arc4random()%(3000))-1500), Vec3D((arc4random()%(2))-1, (arc4random()%(2))-1, (arc4random()%(2))-1), 1);
	 }
	 //*/
	
	
	
	
	
	
	rootNodeBounds = { bodies, numBodies };
	for (int i = 0; i < numBodies; i++)
	{
		bodies[i].bodyKey = computeMortonKey(bodies[i].position, rootNodeBounds.size);
	}
	RadixSortBodies(bodies, numBodies);
	//MergeSortBodiesByMortonKey(numBodies, bodies);
	
	//buildLinearHashedOctreeInPlace(LHTree, bodies, numBodies, rootNodeBounds);
	testBuildLinearHashedOctreeInPlace(LHTree, bodies, numBodies, rootNodeBounds);
	
	
	
	bodiesAccelerations = new Vec3D[numBodies];
	
	
	//
	/*
	 cout<<"\n\n\nrootNodeBounds.size: "<<rootNodeBounds.size;
	 cout<<"\n\n";
	 for (int i = 0; i < numBodies; i++)
	 {
	 cout<<"\n\nBody "<<i;
	 cout<<"\nmass: "<< bodies[i].mass;
	 
	 
	 cout<<"\ndisplayRadius: "<< bodies[i].displayRadius;
	 cout<<"\nposition: ("<< bodies[i].position.x; cout<<", "<< bodies[i].position.y; cout<<", "<< bodies[i].position.z;cout<<" )";
	 
	 cout<<"\nvelocity: ("<< bodies[i].velocity.x; cout<<", "<< bodies[i].velocity.y; cout<<", "<< bodies[i].velocity.z;cout<<" )";
	 
	 }
	 
	 
	 cout<<"\n\n\n\n";
	 //*/
	
	
	panningVelocity = ofVec3f(0, 0, 0);
	zoomVelocity = 0.0f;
	canvas.set(-2 * rootNodeBounds.size, -2 * rootNodeBounds.size, 4 * rootNodeBounds.size, 4 * rootNodeBounds.size);
	maxZoomScale = canvas.width / ofGetWidth();//0.001;
	minZoomScale = ofGetWidth() / canvas.width;//0.00000000000001;
}

//--------------------------------------------------------------
void ofApp::update()
{
	rootNodeBounds = { bodies, numBodies };
	for (int i = 0; i < numBodies; i++)
	{
		bodies[i].bodyKey = computeMortonKey(bodies[i].position, rootNodeBounds.size);
	}
	RadixSortBodies(bodies, numBodies);
	//MergeSortBodiesByMortonKey(numBodies, bodies);
	
	
	
	
	///*
	panningVelocity *= interfaceInertiaDamping;
	zoomVelocity *= interfaceInertiaDamping;
	horizontalRotationVelocity *= interfaceInertiaDamping;
	verticalRotationVelocity *= interfaceInertiaDamping;
	offset += panningVelocity * interfaceInertiaDamping;
	zoomScale += zoomVelocity * interfaceInertiaDamping;
	horizontalAngle += horizontalRotationVelocity;
	verticalAngle += verticalRotationVelocity;
	if (isFocusing)
	{
		// Set offset to target object's position, inversely scaled by current zoom factor
		offset = targetOffset / zoomScale;
		offset += (offset + targetOffset) * focusSpeed;
		if ((targetOffset - offset).length() < 1.0f)  // Check if close enough
		{
			isFocusing = false;
		}
	}
	//*/
}

//--------------------------------------------------------------
void ofApp::draw()
{
	//FrustumPlane planes[6];
	//calculateFrustumPlanes(horizontalAngle, verticalAngle, offset, zoomScale, planes); // Use custom frustum planes based on navigation parameters
	
	//for (size_t i = 0; i < numBodies; i++)
	//{
	//	if (isPointInFrustum(glm::vec3(bodies[i].position.x, bodies[i].position.y, bodies[i].position.z), planes))
	//	{
	//		ofDrawSphere(glm::vec3(bodies[i].position.x, bodies[i].position.y, bodies[i].position.z), bodies[i].mass);
	//	}
	//}
	
	ofSetColor(255);
	
	
	///*
	ofDrawBitmapString("FPS: " + ofToString(ofGetFrameRate(), 2), ofGetWidth() - 200, 45);
	ofDrawBitmapString("MAC: " + ofToString(theta, 2), ofGetWidth() - 200, 60);
	ofDrawBitmapString("numBodies: " + ofToString(numBodies, 2), ofGetWidth() - 200, 75);
	ofDrawBitmapString("Zoom Scale: " + ofToString(zoomScale, 10), ofGetWidth() - 200, 90);
	ofDrawBitmapString("Zoom Level: " + ofToString(currentZoomIndex, 2), ofGetWidth() - 200, 105);
	ofDrawBitmapString("Offset: " + ofToString(offset, 2), ofGetWidth() - 200, 120);
	ofDrawBitmapString("Zoom Increment: " + ofToString(zoomIncrement, 10), ofGetWidth() - 200, 135);
	ofDrawBitmapString("Is focusing: " + ofToString(isFocusing, 2), ofGetWidth() - 200, 150);
	
	//ofDrawBitmapString("Tree Visualize Status: " + ofToString(visualizeTree, 2), ofGetWidth() - 200, 105);
	//*/
	
	
	
	///*
	ofPushMatrix();
	ofTranslate(ofGetWidth() * 0.5, ofGetHeight() * 0.5); // First, translate to center
	ofRotateDeg(ofRadToDeg(horizontalAngle), 0, 1, 0);   // Apply the rotations from the center of the screen // Rotate around Y
	ofRotateDeg(ofRadToDeg(verticalAngle), 1, 0, 0);  // Rotate around X
	ofScale(zoomScale, zoomScale, zoomScale);     // Scale
	
	ofTranslate(-ofGetWidth() * 0.5 + offset.x, -ofGetHeight() * 0.5 + offset.y, offset.z);     // Translate back to original position or offset
																								//*/
	
	
	
	ComputePositionAtHalfTimeStep(dt, bodies, numBodies, bodiesAccelerations);
	
	
	buildLinearHashedOctreeInPlace(LHTree, bodies, numBodies, rootNodeBounds);
	
	if (visualizeTree)
	{
		LHTree.visualizeTree();
	}
	
	
	// For the two HOTNode lists memory is allocated for an array of 'numBodies' pointers to HOTNode object instantiations.
	// However, they do not actually create any HOTNode objects. The array is initialized with pointers,
	// but all of these pointers are uninitialized and do not point to any HOTNode objects yet.
	HOTNode** walkList = new  HOTNode * [numBodies]; // This should be based on the depth of tree, 'DEFAULT_HASHED_OCTREE'
	HOTNode** interactList = new  HOTNode * [numBodies];
	
	//const HOTNode** interactList = (const HOTNode**) malloc(sizeof(HOTNode*)*1000); //  the pointers in the array cannot be modified to point to different HOTNode objects. However, the HOTNode objects themselves can still be modified if they are not declared as const.
	
	
	//Vec3D* bodiesAccelerations;
	//bodiesAccelerations = new Vec3D[numBodies];
	ResetVec3D(bodiesAccelerations, numBodies);
	
	
	
	ComputeHOTOctreeForce(LHTree, bodies, bodiesAccelerations, numBodies, walkList, interactList, theta); //this function computes the accelerations from gravity for all bodies
	ComputeVelocityAndPosition(dt, bodies, numBodies, bodiesAccelerations);
	
	
	
	delete[] walkList;
	delete[] interactList;
	//delete[] bodiesAccelerations;
	LHTree.deleteTree();
	
	
	VisualizeBodies(bodies, numBodies);
	
	ofPopMatrix();
	//*/
}

//--------------------------------------------------------------
void ofApp::exit(){
	
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
	if (key == 'z') //reset zoom
	{
		zoomScale = 1;
		lastZoomScale = zoomScale;
		offset.set(0, 0);
		
		panningVelocity.set(0, 0, 0);  // Reset panning velocity
		horizontalRotationVelocity = 0.0f;
		verticalRotationVelocity = 0.0f;
		zoomVelocity = 0.0f;  // Reset zoom velocity
		
	}
	
	///*
	if (key == 'f')  // Set the target position
	{
		isFocusing = !isFocusing;
		targetOffset = ofVec3f(-ofGetWidth() * 0.5 + ofGetMouseX(), -ofGetHeight() * 0.5 + ofGetMouseY(), 0);
		//zoomScale = 1.0f;
	}
	
	if (key == 'v')
	{
		visualizeTree = !visualizeTree;
	}
	
	
	if(key == 'l') //large view of galaxy//zoomed in on sun
	{
		zoomScale = 0.000000001;
	}
	if(key == 's') //zoomed in on sun
	{
		zoomScale = 0.125;
	}
	
	
	
	if (key == 'o') // Zoom to oort cloud
	{
		zoomScale = 0.0000000001;
		lastZoomScale = zoomScale;
		offset.set(0, 0);
		
		panningVelocity.set(0, 0, 0);  // Reset panning velocity
		horizontalRotationVelocity = 0.0f;
		verticalRotationVelocity = 0.0f;
		zoomVelocity = 0.0f;  // Reset zoom velocity
		
	}
	
	
	
	if (key == OF_KEY_UP) // Scale order of magnitude of zoom by arrow key up
	{
		// Move 'down' in the vector => bigger scale
		currentZoomIndex++;// = (currentZoomIndex <= 14 ? currentZoomIndex++ : currentZoomIndex);
		
		zoomScale = ZOOM_LEVELS[currentZoomIndex];
		zoomIncrement = zoomScale * 0.01;
	}
	
	if (key == OF_KEY_DOWN) // Scale order of magnitude of zoom by arrow key down
	{
		// Move 'up' in the vector => smaller scale
		currentZoomIndex--;// = (currentZoomIndex > 0 ? currentZoomIndex-- : currentZoomIndex);
		
		
		zoomScale = ZOOM_LEVELS[currentZoomIndex];
		zoomIncrement = zoomScale * 0.01;
		
	}
	
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
	
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {
	
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button)
{
	
	
	///*
	ofVec3f mousePos(x, y, 0);  // Z is 0 for basic mouse input
	ofVec3f dragOffset = mousePos - dragStartPos;
	
	
	
	
	
	
	
	
	if (button == 0) //left mouse button for panning
	{
		//
		/* -------- DIRECTIONAL VECTOR PANNING, LESS EXPENSIVE, LESS ACCURATE --------
		 //*--- 2 of 3 three dimensions of panning stop working if the orientation of the scene is changed
		 //*--- updirection has to assume 1 direction is up relative to the other 2, so when this changes it loses it's ability to
		 //*--- calculate panning vectors based on the rotation, the viewer's horizontal and vertical angles have changed
		 ofVec3f rightDirection(cos(horizontalAngle), 0, -sin(horizontalAngle));   // cross product with up vector (0,1,0)
		 ofVec3f upDirection(sin(verticalAngle) * sin(horizontalAngle), cos(verticalAngle), sin(verticalAngle) * cos(horizontalAngle));
		 // New depth direction for forward/backward movement based on viewpoint
		 ofVec3f depthDirection(sin(horizontalAngle), 0, cos(horizontalAngle));
		 // Adjust the offset based on the panning vectors. We're using dragOffset.y for depth.
		 offset += rightDirection * dragOffset.x + upDirection * dragOffset.y + depthDirection * dragOffset.y;
		 //*/
		
		
		
		///* -------- MATRIX INVERSION PANNING, MORE EXPENSIVE, MORE ACCURATE --------
		//*--- Uses a rotation matrix to transform the 2D drag into the 3D world space by applying the inverse of the current view rotation
		
		// Calculate 3D drag offset based on current rotation angles
		ofVec3f dragOffset3D(dragOffset.x, dragOffset.y, 0);
		
		// Create a rotation matrix based on the horizontal and vertical angles
		ofMatrix4x4 rotMatrix;
		rotMatrix.makeIdentityMatrix();
		rotMatrix.rotate(-ofRadToDeg(horizontalAngle), 0, 1, 0);
		rotMatrix.rotate(-ofRadToDeg(verticalAngle), 1, 0, 0);
		
		// Rotate the 3D drag offset by the inverse of the current rotation matrix
		ofVec3f rotatedOffset = dragOffset3D * rotMatrix.getInverse();
		
		// Add the rotated offset to the existing offset
		offset += rotatedOffset;
		//*/
		
	}
	else if (button == 2) //right mouse button for rotating
	{
		// Convert drag movement to spherical coordinates rotation
		double sensitivity = 0.005f;  // Adjust this value for rotation sensitivity
		horizontalRotationVelocity = dragOffset.x * sensitivity;
		verticalRotationVelocity = -dragOffset.y * sensitivity;
		
		horizontalAngle += horizontalRotationVelocity;
		verticalAngle += verticalRotationVelocity;
	}
	
	dragStartPos = mousePos;
}




//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button)
{
	isDragging = true;
	dragStartPos.set(x, y, 0);  // Z is 0 for basic mouse input
	
	
	if (button == 0) //left mouse button for panning
	{
		panningVelocity.set(0, 0, 0);  // Reset panning velocity
		offsetAtMousePress = offset;  // Store the offset value when the mouse is pressed
	}
	else if (button == 2) // right mouse button for rotating
	{
		horizontalRotationVelocity = 0.0;
		verticalRotationVelocity = 0.0;
		
		horizontalAngleAtMousePress = horizontalAngle;
		verticalAngleAtMousePress = verticalAngle;
	}
}



//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button)
{
	isDragging = false;
	if (button == 0)  // Panning
	{
		panningVelocity = (offset - offsetAtMousePress) * 0.1;  // Calculate velocity based on the difference from mouse press to release
	}
	else if (button == 2) // Rotating
	{
		horizontalRotationVelocity = (horizontalAngle - horizontalAngleAtMousePress) * 0.1;  // Calculate rotation velocity based on angle differences
		verticalRotationVelocity = (verticalAngle - verticalAngleAtMousePress) * 0.1;
	}
}







/*
 This function handles mouse scroll events to adjust the zoom level of the view.
 It updates the zoom scale based on the scroll input, ensuring it remains within
 defined limits. The function also adjusts the view offset to zoom towards the
 mouse position, providing a more intuitive zooming experience.
 */
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY)
{
	/// Calculate the amount to adjust the zoom scale based on scroll input
	double scrollAmount = (scrollY) * zoomIncrement;  // Receive amount scrolled, controlled by the arbitrarily defined zoomIncrement
	zoomScale += (scrollAmount);   // Update zoomScale
	
	
	
	///*
	// A: If you want to clamp to the min & max of your *current* level by some margin,
	// you can do so or just skip. For instance, if you want to not exceed 10x difference
	// from ZOOM_LEVELS[currentZoomIndex].
	// e.g.:
	//if(currentZoomIndex > (0) && currentZoomIndex < 14)
	//{
	double baseLevel = ZOOM_LEVELS[currentZoomIndex];
	double lowerLimit = ZOOM_LEVELS[currentZoomIndex+1]; // 10% below
	double upperLimit = ZOOM_LEVELS[currentZoomIndex-1]; // 10x above
	zoomScale = ofClamp(zoomScale, lowerLimit, upperLimit);
	//}
	//else
	//{
	/// Update the maximum and minimum zoom scales based on the canvas and window dimensions
	//	maxZoomScale = (canvas.width*ZOOM_LEVELS[currentZoomIndex]) / ofGetWidth(); // as small as it can get
	//	minZoomScale = ofGetWidth() / (canvas.width*ZOOM_LEVELS[currentZoomIndex]); // as large as it can get
	
	/// Clamp the zoom scale to ensure it remains within valid limits
	//	zoomScale = ofClamp(zoomScale, minZoomScale, maxZoomScale); // Make sure zoomScale is within a valid range
	//}
	//*/
	
	
	
	/// Store the current zoom scale for later comparison
	lastZoomScale = zoomScale;
	zoomVelocity = 0.0;  // Reset zoom velocity
	
	
	
	
	
	
	
	
	/// Adjust the offset to zoom towards the mouse position
	offset -= (ofVec3f(x - offset.x, y - offset.y, offset.z) * (zoomScale - lastZoomScale)) / lastZoomScale; // Zooms to mouse
																											 // Allows user to zoom into the position of the mouse by offsetting the coordinate system relative to the position of the mouse
	
	/// Calculate the zoom velocity for potential use in animations or further adjustments
	zoomVelocity = (zoomScale - lastZoomScale) * 0.1;  // Multiplied by a factor for control
	
	
	// make zoom increment proportional to the current zoom scale to ensure smooth navigation of the scene, where the increment is automatically adjusted even as values of zoom fluctuate using rigorous mathematical functions with their basis in the stability of the rate of zoom, so when zoomed in close the increment will adjust to reflect gradual changes in scope, and when zoomed out the increment will adjust to reflect the broader changes in scope, this will allow for a more intuitive navigation of the scene.
	// Meaning, the rate at which the zoom scale changes will be determined by several factors including the current zoomscale relative to the maximum, minimum, and equilibrium zoomscales
	
	
	// zoom scale will be refactored so that the up and down arrow keys adjust the zoom scale by orders of magnitude of 10(so either 0.1 or 10 times the current zoom scale) to allow for distinct navigational levels enabling both fine and coarse control of the view, while the scroll wheel will then accomodate for more precise navigational inputs within the current zoom level(levels = start from 10x the bounds of the canvas, from here there will be an additional level for every order of magnitude of 10, so 10x, 1x, 0.1x, 0.01x, 0.001x, 0.0001x, 0.00001x, 0.000001x, 0.0000001x, 0.00000001x, 0.000000001x, 0.0000000001x, 0.00000000001x, 0.000000000001x, 0.0000000000001x, 0.00000000000001x, 0.000000000000001x, 0.0000000000000001x, 0.00000000000000001x, 0.000000000000000001x, 0.0000000000000000001x, 0.00000000000000000001x, 0.000000000000000000001x, 0.0000000000000000000001x, 0.00000000000000000000001x, 0.000000000000000000000001x, 0.0000000000000000000000001x, 0.00000000000000000000000001x, 0.000000000000000000000000001x, 0.0000000000000000000000000001x, and so on until the minimum zoom scale is reached)
	// the up and down arrow keys will adjust the zoom scale by orders of magnitude of 10 to allow for distinct navigational levels enabling both fine and coarse control of the view, while the scroll wheel will adjust the zoom scale by a factor of 1.1 or 0.9, relative to the current level, to allow for dynamic navigation via scroll wheel within levels.
	
	
	
	
	
	
	
}






//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){
	
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){
	
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h)
{
	maxZoomScale = canvas.width / ofGetWidth();
	minZoomScale = ofGetWidth() / canvas.width;
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
	
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){
	
}



































