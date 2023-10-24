#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <functional>

#include "Utils/Globals.h"  // #include <glad/glad.h>

#include "../util.hpp"
#include "Utils/Utils.h"

namespace gfx {
static std::string shader_prepend = "#version 430 core\n";
static std::string shader_root = "src/sim-fluid-loganzartman/shader";

class Program {
public:
    Program() : name("[unnamed]") {}
    Program(std::string name) : name(name) {}
    ~Program() {
        if (id) { GlobalsGL::f()->glDeleteProgram(id); }
        if (vertex_id) { GlobalsGL::f()->glDeleteShader(vertex_id); }
        if (geometry_id) { GlobalsGL::f()->glDeleteShader(geometry_id); }
        if (fragment_id) { GlobalsGL::f()->glDeleteShader(fragment_id); }
        if (compute_id) { GlobalsGL::f()->glDeleteShader(compute_id); }
    }

    Program& vertex(std::initializer_list<std::string> srcs) {
        compile_shader(srcs, vertex_id, GL_VERTEX_SHADER, "vertex");
        return *this;
    }

    Program& geometry(std::initializer_list<std::string> srcs) {
        compile_shader(srcs, geometry_id, GL_GEOMETRY_SHADER, "geometry");
        return *this;
    }

    Program& fragment(std::initializer_list<std::string> srcs) {
        compile_shader(srcs, fragment_id, GL_FRAGMENT_SHADER, "fragment");
        return *this;
    }

    Program& compute(std::initializer_list<std::string> srcs) {
        compile_shader(srcs, compute_id, GL_COMPUTE_SHADER, "compute");
        return *this;
    }

    void compile_shader(std::initializer_list<std::string> srcs, GLuint& dest, GLenum type, std::string type_name) {
        std::string src = read_sources(srcs);
        dest = GlobalsGL::f()->glCreateShader(type);
        const char* src_str = src.c_str();
        GlobalsGL::f()->glShaderSource(dest, 1, &src_str, NULL);
        GlobalsGL::f()->glCompileShader(dest);
        try {
            check_shader_errors(dest);
        } catch (std::runtime_error& e) {
            throw std::runtime_error(type_name + " shader compilation error in files: " + join(srcs.begin(), srcs.end(), ", ") + "\n" + e.what());
        }
    }

    Program& compile() {
        id = GlobalsGL::f()->glCreateProgram();

        // attach shaders
        if (compute_id) {
            GlobalsGL::f()->glAttachShader(id, compute_id);
        }
        else {
            if (!vertex_id) { throw std::runtime_error("Compiling program without vertex shader loaded."); }
            if (!fragment_id) { throw std::runtime_error("Compiling program without fragment shader loaded."); }
            GlobalsGL::f()->glAttachShader(id, vertex_id);
            if (geometry_id) { GlobalsGL::f()->glAttachShader(id, geometry_id); }
            GlobalsGL::f()->glAttachShader(id, fragment_id);
        }
        
        GlobalsGL::f()->glLinkProgram(id);
        check_program_errors(id);
        return *this;
    }

    GLint uniform_loc(std::string uname) {
        GLint location = GlobalsGL::f()->glGetUniformLocation(id, uname.c_str());
        // if (location < 0) { std::cerr << "Warning in " << name << ": invalid or unused uniform: " << uname << std::endl; }
        return location;
    }
    
    void use() {
        if (!id) { throw std::runtime_error("Trying to use program that is not compiled."); }
        GlobalsGL::f()->glUseProgram(id);
    }

    void validate() {
        GlobalsGL::f()->glValidateProgram(id);
        check_validation(id);
    }

    void disuse() {
        GlobalsGL::f()->glUseProgram(0);
    }

    std::string name;
    GLuint id = 0;
    GLuint vertex_id = 0;
    GLuint geometry_id = 0;
    GLuint fragment_id = 0;
    GLuint compute_id = 0;

private:
    /**
     * Read and concatenate shader sources from the shader_root
     */
    std::string read_sources(std::initializer_list<std::string> srcs) {
        std::stringstream result;
        result << shader_prepend;
        for (auto& src : srcs) {
            result << file_read(shader_root + "/" + src) << std::endl;
        }
        return result.str();
    }

    void check_validation(GLuint shader) {
        GLint is_ok = 0;
        GlobalsGL::f()->glGetProgramiv(id, GL_VALIDATE_STATUS, &is_ok);
        if (!is_ok) {
			GLint max_length = 0;
            GlobalsGL::f()->glGetProgramiv(id, GL_INFO_LOG_LENGTH, &max_length);

            // The maxLength includes the NULL character
            std::vector<GLchar> info_log(max_length);
            GlobalsGL::f()->glGetProgramInfoLog(id, max_length, NULL, &info_log[0]);

            const std::string err = "GLSL validation error for program: ";
            throw std::runtime_error(err + name + "\n" + info_log.data());
		}
    }

    void check_shader_errors(GLuint shader) {
        GLint is_ok = 0;
        GlobalsGL::f()->glGetShaderiv(shader, GL_COMPILE_STATUS, &is_ok);
        if (!is_ok) {
			GLint max_length = 0;
            GlobalsGL::f()->glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &max_length);

            // The maxLength includes the NULL character
            std::vector<GLchar> info_log(max_length);
            GlobalsGL::f()->glGetShaderInfoLog(shader, max_length, NULL, &info_log[0]);

            const std::string err = "GL Program Validation error: \n";
            throw std::runtime_error(err + info_log.data());
		}
    }

    void check_program_errors(GLuint shader) {
        GLint is_ok = 0;
        GlobalsGL::f()->glGetProgramiv(id, GL_LINK_STATUS, &is_ok);
        if (!is_ok) {
			GLint max_length = 0;
            GlobalsGL::f()->glGetProgramiv(id, GL_INFO_LOG_LENGTH, &max_length);

            // The maxLength includes the NULL character
            std::vector<GLchar> info_log(max_length);
            GlobalsGL::f()->glGetProgramInfoLog(id, max_length, NULL, &info_log[0]);

            const std::string err = "GLSL linking error for program: ";
            throw std::runtime_error(err + name + "\n" + info_log.data());
		}
    }
};
}
