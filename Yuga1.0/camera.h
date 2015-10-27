#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>

#include <iostream>
using namespace std;


#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

class CameraBase
{
public:
	CameraBase();
	CameraBase(float R, float alpha, float theta, float wWidth, float wHeight);

	void                SetView(glm::vec3 pEye, glm::vec3 pLookat);
	void                SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane);


	glm::mat4 getMVP()                                    { return mProj * mView * mWorld; }
	glm::mat4 getM()                                      { return mWorld; }
	glm::mat4 getV()                                      { return mView; }

	void initMousePosition(float x, float y)              { SetCurMousePosition(x, y); SetPreMousePosition(x, y); }
	void SetMouseLButtonStat(bool stat)                   { mbMouseLButtonDown = stat; }
	void SetMouseLWheelStat(bool stat)                    { mbMouseWheelRoll = stat; }
	void SetMouseRButtonStat(bool stat)                   { mbMouseRButtonDown = stat; }

	void SetCurMousePosition(float x, float y)            { curMousePosition.x = x; curMousePosition.y = y; }
	void SetPreMousePosition(float x, float y)            { preMousePosition.x = x; preMousePosition.y = y; }
	glm::vec2 GetCurMousePosition()            { return curMousePosition; }
	glm::vec2 GetPreMousePosition()            { return preMousePosition; }

	bool IsMouseLButtonDown() const                       { return mbMouseLButtonDown; }
	bool IsMouseWheelRoll() const                         { return mbMouseRButtonDown; }
	bool IsMouseRButtonDown() const                       { return mbMouseRButtonDown; }

	void ComputeCameraPosition();
	virtual void MoveRight(double amount);
	virtual void MoveUp(double amount);
	virtual void ZoomIn(double amount);

	void Reset();

protected:
	glm::vec3 mFocus;
	glm::vec3 mCameraPosition;
	glm::vec3 mUp;

	float     mRadius;
	float     mAlpha;
	float     mTheta;

	float     movementSensitivity;

	glm::mat4 mView;
	glm::mat4 mProj;
	glm::mat4 mWorld;

	bool mbMouseLButtonDown;    // True if left button is down 
	bool mbMouseWheelRoll;          // True if middle wheel is roll 
	bool mbMouseRButtonDown;    // True if right button is down 

	glm::vec2 curMousePosition;
	glm::vec2 preMousePosition;
};

class ObjectCamera : public CameraBase
{
public:
	ObjectCamera(glm::vec3 cPosition, glm::vec3 fPosition, float wWidth, float wHeight);
	~ObjectCamera();

	void MoveRight(double amount);
	void MoveUp(double amount);
	void ZoomIn(double amount);
	

	glm::vec3                                             ConvertXY(float x, float y);
	void SetRotation();

	void setMmworldScle()                                 { mWorld = mScale * mWorld; };

	void setScaleFactor(float x)                          { mScale = glm::scale(glm::mat4(), glm::vec3(x)); }
	void setWH(float w, float h)                          { mWidth = w, mHeight = h; }

	void computeAxis()
	{
		zAxis = mCameraPosition - mFocus;
		zAxis = glm::normalize(zAxis);
		xAxis = glm::cross(zAxis, mUp);
		xAxis = glm::normalize(xAxis);
		yAxis = glm::cross(zAxis, xAxis);
		yAxis = glm::normalize(yAxis);
	}

	void camara2worldVector(double x, double y, double z, double *w);

private:
	glm::mat4 mScale;

	float mWidth;
	float mHeight;

	glm::quat RotationUp;//rotate quaternion
	glm::quat RotationRight;

	glm::vec3 xAxis, yAxis, zAxis;
};

#endif 