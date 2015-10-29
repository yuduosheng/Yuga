#include "App.h"
#include "camera.h"
#include "tetMesh.h"
#include "FULLIntegrator.h"
#include "ProjectAndUnProject.h"
#include "TIMING_BREAKDOWN.h"
namespace shader
{
	GLuint loadShader(const char * filename, GLenum shader_type, bool check_errors)
	{
		GLuint result = 0;
		string shaderCode = "";
		ifstream codeStream(filename, ios::in);

		if (codeStream.is_open())
		{
			string Line = "";
			while (getline(codeStream, Line))
			{
				shaderCode += Line;
				shaderCode += "\n";
			}

			codeStream.close();
		}
		else
		{
			cout << filename << " can not open! " << endl;
		}
		result = glCreateShader(shader_type);

		// Compile Vertex Shader
		//cout << "Compiling shader : " << filename << endl;
		char const * codePointer = shaderCode.c_str();

		glShaderSource(result, 1, &codePointer, NULL);
		glCompileShader(result);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(result, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(result, 4096, NULL, buffer);

				cout << filename << buffer << endl;
			}
		}

		return result;
	}

}


class Test : public App
{
public:
	Test();
	~Test();

	bool                    Init();
	void                    UpdateScene();
	void                    Rendering();
	void                    onResize(GLFWwindow* window, int w, int h);

	void                    onMouseWheel(GLFWwindow* window, double x, double y);
	void                    onMouseMove(GLFWwindow* window, double xd, double yd);
	void                    onMouseButton(GLFWwindow* window, int button, int action, int mods);
	void                    onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

private:
	void                    buildGeometryBuffers();
	void                    buildShader();
	void                    setGridCellSize();

private:
	ObjectCamera            camera;
	GLuint                  vao;
	GLuint                  program;
	// This will identify our vertex buffer
	GLuint                  vertexbuffer;
	GLuint                  elementbuffer;

	GLuint                  mvp_matrix;
	GLuint                  m_matrix;
	GLuint                  v_matrix;
	GLuint                  l_position;

	TetMesh                 *modelTetMesh;
	FullIntegrator          *integrator;
	bool					isFill;

	int dragStartX, dragStartY;
	int pulledVertex;

	VEC3F vertexForce;
};

int main(void)
{
	Test *theApp = new Test;

	if (!theApp->Init())
		return 0;
	theApp->Run();
	delete theApp;
	return 0;
}

Test::Test() : App(), camera(glm::vec3(0.0f, 2.0f, 10.0f), glm::vec3(0.0f, 2.0f, 0.0f), mWidth, mHeight),
isFill(true)
{
	pulledVertex = -1;
}

Test::~Test()
{
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	// Cleanup VBO and shader
	glDeleteProgram(program);
	glDeleteVertexArrays(1, &vao);
	if (modelTetMesh)
		delete modelTetMesh;
	if (integrator)
		delete integrator;

}

bool Test::Init()
{
	if (!App::Init())
		return false;


	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	buildShader();
	//buildGeometryBuffers();
	modelTetMesh = new TetMesh("capsule2.1");
	integrator = new FullIntegrator(modelTetMesh);

	return true;
}

void Test::onResize(GLFWwindow* window, int nw, int nh)
{
	App::onResize(window, nw, nh);

	glViewport(0, 0, nw, nh);
	camera.setWH(nw, nh);
	camera.SetProj(45.0f, (GLfloat)nw / (GLfloat)nh, 0.01f, 100.0f);

}

void Test::UpdateScene()
{
	TIMING_BREAKDOWN::startFrame();
	static const GLfloat background[] = { 0.1f, 0.1f, 0.1f, 0.1f };
	static const GLfloat one = 1.0f;

	glViewport(0, 0, mWidth, mHeight);
	glClearBufferfv(GL_COLOR, 0, background);
	glClearBufferfv(GL_DEPTH, 0, &one);

	glUseProgram(program);

	glm::mat4 MVP = camera.getMVP();
	glm::mat4 M = camera.getM();
	glm::mat4 V = camera.getV();

	glUniformMatrix4fv(mvp_matrix, 1, GL_FALSE, glm::value_ptr(MVP));
	glUniformMatrix4fv(v_matrix, 1, GL_FALSE, glm::value_ptr(V));
	glUniformMatrix4fv(m_matrix, 1, GL_FALSE, glm::value_ptr(M));
	glm::vec3 lightPos = glm::vec3(10.0f, 10.0f, 10.0f);
	glUniform3f(l_position, lightPos.x, lightPos.y, lightPos.z);

	if (isFill)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	if (pulledVertex != -1)
	{
		integrator->resetExternalForce();
		integrator->addVertexForce(pulledVertex, vertexForce);
	}

	integrator->DoTimeStep();
	TIMING_BREAKDOWN::endFrame();
}
void Test::Rendering()
{
	modelTetMesh->RenderModel();

	/* Draw a triangle 
	glBegin(GL_TRIANGLES);

	glColor3f(1.0, 0.0, 0.0);    // Red
	glVertex3f(0.0, 1.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);    // Green
	glVertex3f(-1.0, -1.0, 0.0);

	glColor3f(0.0, 0.0, 1.0);    // Blue
	glVertex3f(1.0, -1.0, 0.0);

	glEnd();*/
	// 1rst attribute buffer : vertices
	/*
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
		);

	// Draw the triangle !
	//glDrawArrays(GL_TRIANGLES, 0, 3); // Starting from vertex 0; 3 vertices total -> 1 triangle

	// Index buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

	// Draw the triangles !
	glDrawElements(
		GL_TRIANGLES,      // mode
		12,    // count
		GL_UNSIGNED_SHORT,   // type
		(void*)0           // element array buffer offset
		);
*/
}
void Test::buildGeometryBuffers()
{
	/*
	// An array of 3 vectors which represents 3 vertices
	static const GLfloat g_vertex_buffer_data[] = {
		-0.70710678,  -0.40824829,  0.57735027,
		0.00000000,   0.81649658,   0.57735027,
		0.70710678,   -0.40824829,  0.57735027,
		-0.00000000,  -0.81649658,  -0.57735027,
		-0.70710678,  0.40824829,   -0.57735027,
		0.70710678,   0.40824829,   -0.57735027
	};
	static const GLushort g_index_buffer_data[] = {
		1, 2, 5,
		1, 5, 4,
		1, 4, 3,
		1, 3, 2,
		2, 3, 6,
		2, 6, 5,
		3, 4, 6,
		4, 5, 6
	};*/
	// An array of 3 vectors which represents 3 vertices
	static const GLfloat g_vertex_buffer_data[] = {
		0, 0, 0,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};
	static const GLushort g_index_buffer_data[] = {
		1, 0, 2,
		0, 3, 2,
		1, 0, 3,
		1, 2, 3
	};
	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);

	// The following commands will talk about our 'vertexbuffer' buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

	// Give our vertices to OpenGL.
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

	// Generate a buffer for the indices

	glGenBuffers(1, &elementbuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(g_index_buffer_data) * sizeof(unsigned short), g_index_buffer_data, GL_STATIC_DRAW);

}
void Test::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::loadShader("lighting.vs.glsl", GL_VERTEX_SHADER, true);
	fs = shader::loadShader("lighting.fs.glsl", GL_FRAGMENT_SHADER, true);

	//vs = shader::loadShader("ColorVertexShader.glsl", GL_VERTEX_SHADER, true);
	//fs = shader::loadShader("ColorFragmentShader.glsl", GL_FRAGMENT_SHADER, true);
	if (program)
		glDeleteProgram(program);

	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);

	glLinkProgram(program);

	glDeleteShader(vs);
	glDeleteShader(fs);

	mvp_matrix = glGetUniformLocation(program, "MVP");
	v_matrix = glGetUniformLocation(program, "V");
	m_matrix = glGetUniformLocation(program, "M");
	l_position = glGetUniformLocation(program, "LightPosition_worldspace");

}

void Test::onMouseWheel(GLFWwindow* window, double x, double y)
{
	float scale = 1.0f;
	float mouseWheelScale = 0.1f;
	scale += mouseWheelScale  * (float)y;
	camera.setScaleFactor(scale);
	camera.setMmworldScle();
}
void Test::onMouseMove(GLFWwindow* window, double xd, double yd)
{
	double x = xd;
	double y = yd;
	if (camera.IsMouseLButtonDown())
	{
		camera.SetCurMousePosition(x, y);
		camera.SetRotation();
		camera.SetPreMousePosition(x, y);

		if (pulledVertex != -1)
		{
			double forceX = (x - dragStartX);
			double forceY = -(y - dragStartY);

			double externalForce[3];
			camera.camara2worldVector(forceX, forceY, 0, externalForce);

			vertexForce[0] = externalForce[0];
			vertexForce[1] = externalForce[1];
			vertexForce[2] = externalForce[2];
		}
	}

}
void Test::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;

	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{
		GLdouble model[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, model);

		GLdouble proj[16];
		glGetDoublev(GL_PROJECTION_MATRIX, proj);

		GLint view[4];
		glGetIntegerv(GL_VIEWPORT, view);

		glfwGetCursorPos(window, &xd, &yd);

		int width, height;
		glfwGetFramebufferSize(window, &width, &height);

		int winX = xd;
		int winY = height - yd;

		float zValue;
		glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zValue);

		GLubyte stencilValue;
		glReadPixels(winX, winY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, &stencilValue);
		if (stencilValue == 1)
		{
			GLdouble worldCoor[3];
			_glUnProject(winX, winY, zValue, model, proj, view, worldCoor);

			dragStartX = xd;
			dragStartY = yd;
			VEC3F pos(worldCoor[0], worldCoor[1], worldCoor[2]);

			pulledVertex = modelTetMesh->getClosestVertex(pos);

			printf("Clicked on vertex: %d (0-indexed)\n", pulledVertex);
		}
		else
		{
			printf("Clicked on empty stencil: %d.\n", stencilValue);
		}
		if (pulledVertex == -1)
		{
			camera.SetMouseLButtonStat(true);
			glfwGetCursorPos(window, &xd, &yd);
			camera.initMousePosition(xd, yd);
		}
	}
	else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE))
	{
		camera.SetMouseLButtonStat(false);
		pulledVertex = -1;
	}

}
void Test::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	if ((key == GLFW_KEY_1) && (action == GLFW_PRESS))
	{
		isFill = !isFill;
	}
	if ((key == GLFW_KEY_T) && (action == GLFW_PRESS))
	{
		TIMING_BREAKDOWN::printTimingBreakdown();
	}
}