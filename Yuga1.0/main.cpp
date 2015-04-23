#include "App.h"



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

	GLuint                  vao;
	GLuint                  program;
	// This will identify our vertex buffer
	GLuint                  vertexbuffer;
	//GLuint                  mvp_matrix;
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

Test::Test() : App()
{

}

Test::~Test()
{
	// Cleanup VBO and shader
}

bool Test::Init()
{
	if (!App::Init())
		return false;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	buildShader();
	buildGeometryBuffers();
	return true;
}

void Test::onResize(GLFWwindow* window, int nw, int nh)
{
	App::onResize(window, nw, nh);

	glViewport(0, 0, nw, nh);

}

void Test::UpdateScene()
{
}
void Test::Rendering()
{
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
	glDrawArrays(GL_TRIANGLES, 0, 3); // Starting from vertex 0; 3 vertices total -> 1 triangle

}
void Test::buildGeometryBuffers()
{
	// An array of 3 vectors which represents 3 vertices
	static const GLfloat g_vertex_buffer_data[] = {
		-1.0f, -1.0f, 0.0f,
		1.0f, -1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
	};

	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);

	// The following commands will talk about our 'vertexbuffer' buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

	// Give our vertices to OpenGL.
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

}
void Test::buildShader()
{
	GLuint vs;
	GLuint fs;

	vs = shader::loadShader("ColorVertexShader.glsl", GL_VERTEX_SHADER, true);
	fs = shader::loadShader("ColorFragmentShader.glsl", GL_FRAGMENT_SHADER, true);

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

	//mvp_matrix = glGetUniformLocation(program, "MVP");

}

void Test::onMouseWheel(GLFWwindow* window, double x, double y)
{

}
void Test::onMouseMove(GLFWwindow* window, double xd, double yd)
{
	double x = xd;
	double y = yd;
}
void Test::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;

	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{

	}

}
void Test::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS))
	{
	}

}