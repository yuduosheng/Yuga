#include "camera.h"
CameraBase::CameraBase()
{

}
CameraBase::CameraBase(float R, float alpha, float theta, float wWidth, float wHeight)
{
	mFocus = glm::vec3(0.0f, 0.0f, 0.0f);
	mUp = glm::vec3(0.0f, 1.0f, 0.0f);
	movementSensitivity = 0.01f;
	mWorld = glm::mat4(1.0f);

	mRadius = R;
	mAlpha = alpha;
	mTheta = theta;

	ComputeCameraPosition();
	SetProj(45.0f, (float)wWidth / wHeight, 0.1f, 100.f);

	mbMouseLButtonDown = false;
	mbMouseWheelRoll = false;
	mbMouseRButtonDown = false;
}

void CameraBase::ComputeCameraPosition()
{
	//cameraPosition = focusPosition + 
	//  Vec3(R * cos(Phi) * cos (Theta), R * sin(Theta), -R * sin(Phi) * cos (Theta));
		mCameraPosition = glm::vec3(
			mFocus.x + mRadius * cos(mAlpha) * cos(mTheta),
			mFocus.y + mRadius * sin(mTheta),
			mFocus.z + mRadius * sin(mAlpha) * cos(mTheta)
			);


	SetView(mCameraPosition, mFocus);
}

void CameraBase::SetView(glm::vec3 eye, glm::vec3 target)
{
	mView = glm::lookAt(eye, target, mUp);
}

void CameraBase::SetProj(float fFov, float fAspect, float fNearPlane, float fFarPlane)
{
	mProj = glm::perspective(fFov, fAspect, fNearPlane, fFarPlane);
}

void CameraBase::MoveRight(double amount)
{
	mAlpha += amount * movementSensitivity;
	ComputeCameraPosition();
}

void CameraBase::MoveUp(double amount)
{
	mTheta += amount * movementSensitivity;
	cout << mTheta << endl;
	if (mTheta > 89.0 * M_PI / 180)
		mTheta = 89.0 * M_PI / 180;

	if (mTheta < -89.0 * M_PI / 180)
		mTheta = -89.0 * M_PI / 180;

	ComputeCameraPosition();

}

void CameraBase::ZoomIn(double amount)
{
	mRadius -= amount * movementSensitivity;

	/*
	if (R < fabs(movementSensitivity))
	R = fabs(movementSensitivity);
	*/

	ComputeCameraPosition();

}


void CameraBase::Reset() // resets position to defaults
{
	mRadius = 2.0;
	mAlpha = 0.0;
	mTheta = 0.0;
	ComputeCameraPosition();
}


ObjectCamera::ObjectCamera(glm::vec3 cPosition, glm::vec3 fPosition, float wWidth, float wHeight)
{
	mUp = glm::vec3(0.0f, 1.0f, 0.0f);

	mWidth = wWidth;
	mHeight = wHeight;
	mRadius = glm::min(wWidth / 2, wHeight / 2);

	movementSensitivity = 0.1f;
	mWorld = glm::mat4(1.0f);
	


	mCameraPosition = cPosition;
	mFocus = fPosition;

	RotationRight = glm::quat();
	RotationUp = glm::quat();


	mView = glm::lookAt(mCameraPosition, mFocus, mUp);
	SetProj(45.0f, (float)wWidth / wHeight, 0.1f, 100.f);

	mbMouseLButtonDown = false;
	mbMouseWheelRoll = false;
	mbMouseRButtonDown = false;
}

ObjectCamera::~ObjectCamera()
{
}
void ObjectCamera::MoveRight(double amount)
{
	glm::vec3 c2f = mCameraPosition - mFocus;
	glm::vec3 x = glm::cross(c2f, mUp);
	glm::vec3 y = glm::cross(x, c2f);

	glm::normalize(y);

	float angle = amount * movementSensitivity * (M_PI / 180);

	RotationRight = glm::angleAxis(angle, y);
}
void ObjectCamera::MoveUp(double amount)
{
	glm::vec3 c2f = mCameraPosition - mFocus;
	glm::vec3 x = glm::cross(c2f, mUp);

	glm::normalize(x);

	float angle = amount * movementSensitivity * (M_PI / 180);

	RotationUp = glm::angleAxis(-angle, x);
}
void ObjectCamera::ZoomIn(double amount)
{
}
glm::vec3 ObjectCamera::ConvertXY(float xx, float yy)
{
	float x = xx - mWidth / 2;
	float y = mHeight / 2 - yy;

	double d = x*x + y*y;
	double radiusSquared = mRadius * mRadius;
	if (d > radiusSquared)
	{
		return glm::vec3((float)x, (float)y, 0);
	}
	else
	{
		return glm::vec3((float)x, (float)y, sqrt(radiusSquared - d));
	}
}
void ObjectCamera::SetRotation()
{
	glm::vec3 preV, curV, axis;

	preV = ConvertXY(preMousePosition.x, preMousePosition.y);
	curV = ConvertXY(curMousePosition.x, curMousePosition.y);

	axis = glm::cross(preV, curV);
	axis = glm::normalize(axis);

	float val = glm::dot(preV, curV);
	val = val / (glm::length(preV) * glm::length(curV))-0.00001;
	val = acos(val);

	glm::mat4 rotateM = glm::rotate(val*15, axis);

	mWorld = rotateM * mWorld;
}