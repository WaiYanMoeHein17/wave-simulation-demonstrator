#ifndef VISUALISATION_HPP
#define VISUALISATION_HPP

#include "Shader.hpp"
#include <glm/glm.hpp>

#include "util.hpp"
#include "ColorMap.hpp"

class Visualisation {
private:
public:
    cv::Mat terrainImage;
    GLuint terrainTexture;

    cv::Mat terrainSavedImage; // for remote

    cv::Mat waterImage;
    GLuint waterTexture;
    bool useWaterTexture = false;

    int colorMapIndex = 0;

    GLuint VAO, VBO, shaderProgram;
    float quadVertices[24] = {
        // Positions     // Texture Coords
        -1.0f,  1.0f,   0.0f, 1.0f,  // Top-left
        -1.0f, -1.0f,   0.0f, 0.0f,  // Bottom-left
         1.0f, -1.0f,   1.0f, 0.0f,  // Bottom-right
    
        -1.0f,  1.0f,   0.0f, 1.0f,  // Top-left
         1.0f, -1.0f,   1.0f, 0.0f,  // Bottom-right
         1.0f,  1.0f,   1.0f, 1.0f   // Top-right
    };

    float contourLineFactor;
    bool useGradientColor;
    bool useGrayscale;
    bool paused;

    float simulationScale;
    float simulationOffset;

    Visualisation() : terrainTexture(-1), VAO(-1), VBO(-1), shaderProgram(-1),
                      contourLineFactor(14.0f), useGradientColor(true), useGrayscale(false),
                      simulationScale(1.0f), simulationOffset(0.0f), paused(false) {

        // compile shaders
        unsigned int vertexShader = compileShaderFromFile(GL_VERTEX_SHADER, "../shaders/terrain.vs");
        unsigned int fragmentShader = compileShaderFromFile(GL_FRAGMENT_SHADER, "../shaders/terrain.fs");
        
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);


        // Flip quadVertices vertically
        for (int i = 1; i < 24; i += 4) {
            quadVertices[i] = -quadVertices[i];
        }

        // Create VAO and VBO
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        
        glBindVertexArray(VAO);
        
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
        
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }
    ~Visualisation() {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
    }
    void togglePause() {
        paused = !paused;
    }
    bool isPaused() {
        return paused;
    }
    void toggleGradientColor() {
        useGradientColor = !useGradientColor;
    }
    void toggleGrayscale() {
        useGrayscale = !useGrayscale;
    }
    void incrementContourLineFactor(float increment = 1.0f) {
        contourLineFactor += increment;
    }
    void resetContourLineFactor() {
        contourLineFactor = 14.0f;
    }
    void setContourLineFactor(float factor) {
        contourLineFactor = factor;
    }
    void decrementContourLineFactor(float decrement = 1.0f) {
        contourLineFactor -= decrement;
        if (contourLineFactor < 0.0f) {
            contourLineFactor = 0.0f;
        }
    }
    void toggleWaterTexture() {
        useWaterTexture = !useWaterTexture;
    }

    void draw(GLuint simulationHeightMap = 0, cv::Mat *waterMap = nullptr, bool showSimulation = false) {
        // GLenum error = glGetError();
        // if (error != GL_NO_ERROR) {
        //     std::cerr << "OpenGL Error: " << error << std::endl;
        // }

        if (waterMap != nullptr) {
            cvMat2GLtexture(*waterMap, waterTexture);
        }

        glUseProgram(shaderProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, terrainTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "terrain"), 0);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, simulationHeightMap);
        glUniform1i(glGetUniformLocation(shaderProgram, "waterHeightMap"), 1);

        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, waterTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "waterJetMap"), 2);

        // glm::vec3 geoColorMap[8] = {
        //     glm::vec3(0.0f, 0.0f, 1.0f),  // Deep blue
        //     glm::vec3(0.0f, 0.5f, 1.0f),  // Light blue
        //     glm::vec3(0.0f, 1.0f, 0.0f),  // Green
        //     glm::vec3(1.0f, 1.0f, 0.0f),  // Yellow
        //     glm::vec3(1.0f, 0.5f, 0.0f),  // Orange
        //     glm::vec3(1.0f, 0.0f, 0.0f),  // Red
        //     glm::vec3(0.5f, 0.0f, 0.0f),  // Dark red
        //     glm::vec3(0.5f, 0.5f, 0.5f)   // Gray
        // };
        // glUniform3fv(glGetUniformLocation(shaderProgram, "colorMap"), 8, &geoColorMap[0][0]);
        glUniform3fv(glGetUniformLocation(shaderProgram, "colorMap"), 8, &COLOR_MAPS[colorMapIndex][0][0]);

        glUniform1f(glGetUniformLocation(shaderProgram, "contourLineFactor"), contourLineFactor);
        glUniform1i(glGetUniformLocation(shaderProgram, "gradientColor"), useGradientColor);

        glUniform1i(glGetUniformLocation(shaderProgram, "grayscale"), useGrayscale);

        glUniform1f(glGetUniformLocation(shaderProgram, "simulationScale"), simulationScale);
        glUniform1f(glGetUniformLocation(shaderProgram, "simulationOffset"), simulationOffset);

        glUniform1i(glGetUniformLocation(shaderProgram, "useWaterTexture"), useWaterTexture);
        glUniform1i(glGetUniformLocation(shaderProgram, "showSimulation"), showSimulation);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }

    void updateQuadVertices() {
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(quadVertices), quadVertices);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        for (int i = 0; i < 6; ++i) {
            std::cout << "Vertex " << i << ": (" << quadVertices[i * 4] << ", " << quadVertices[i * 4 + 1] << ")" << std::endl;
        }
        std::cout << std::endl;
    }
    void flipX() {
        for (int i = 0; i < 6; ++i) {
            quadVertices[i * 4] = -quadVertices[i * 4];
        }
        updateQuadVertices();
    }
    void flipY() {
        for (int i = 0; i < 6; ++i) {
            quadVertices[i * 4 + 1] = -quadVertices[i * 4 + 1];
        }
        updateQuadVertices();
    }
    void updateVertexFromMouseInput(GLFWwindow* window, int key, int action) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        int width, height;
        glfwGetWindowSize(window, &width, &height);

        float ndc_x = (xpos / width) * 2.0f - 1.0f;
        float ndc_y = 1.0f - (ypos / height) * 2.0f; // Invert Y

        switch (key) {
            case GLFW_KEY_1:
                quadVertices[0] = ndc_x;
                quadVertices[1] = ndc_y;
                quadVertices[12] = ndc_x;
                quadVertices[13] = ndc_y;
                break;
            case GLFW_KEY_3:
                quadVertices[4] = ndc_x;
                quadVertices[5] = ndc_y;
                break;
            case GLFW_KEY_4:
                quadVertices[8] = ndc_x;
                quadVertices[9] = ndc_y;
                quadVertices[16] = ndc_x;
                quadVertices[17] = ndc_y;
                break;
            case GLFW_KEY_2:
                quadVertices[20] = ndc_x;
                quadVertices[21] = ndc_y;
                break;
        }

        updateQuadVertices();
    }

    void terrainToGL(const cv::Mat& terrain) {
        if (paused) return;
        terrainImage = terrain.clone();
        if (terrainTexture == -1) {
            glGenTextures(1, &terrainTexture);
            std::cout << "Generated texture ID: " << terrainTexture << std::endl;
            std::cout << "Terrain image size: " << terrain.cols << "x" << terrain.rows 
                      << ", Channels: " << terrain.channels() 
                      << ", Depth: " << terrain.depth() << std::endl;
        }
        glBindTexture(GL_TEXTURE_2D, terrainTexture);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        if (terrain.channels() == 3) {
            // RGB image
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, terrain.cols, terrain.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, terrain.data);
        } else if (terrain.channels() == 1) {
            // Grayscale image
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, terrain.cols, terrain.rows, 0, GL_RED, GL_UNSIGNED_BYTE, terrain.data);
        } else if (terrain.depth() == CV_16U) {
            // Z16 image (16-bit depth)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R16, terrain.cols, terrain.rows, 0, GL_RED, GL_UNSIGNED_SHORT, terrain.data);
        } else {
            std::cerr << "Unsupported terrain format: channels=" << terrain.channels() 
                  << ", depth=" << terrain.depth() << std::endl;
        }
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, blurred_grayscale_image.cols, blurred_grayscale_image.rows, 0, GL_RED, GL_UNSIGNED_BYTE, blurred_grayscale_image.data);
        glGenerateMipmap(GL_TEXTURE_2D);

        // Unbind the texture
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void saveImage(const std::string& filename) {
        int width, height;
        glfwGetFramebufferSize(glfwGetCurrentContext(), &width, &height);
        unsigned char* pixels = new unsigned char[width * height * 4];
        glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

        cv::Mat image(height, width, CV_8UC4, pixels);
        cv::Mat flippedImage;
        cv::flip(image, flippedImage, 0);
        cv::Mat rgbImage;
        cv::cvtColor(flippedImage, rgbImage, cv::COLOR_RGBA2BGR);
        cv::imwrite(filename, rgbImage);
        delete[] pixels;

        std::cout << "Image saved as " << filename << std::endl;
    }

    cv::Mat getGLImage() {
        int width, height;
        glfwGetFramebufferSize(glfwGetCurrentContext(), &width, &height);
        unsigned char* pixels = new unsigned char[width * height * 4];
        glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

        cv::Mat image(height, width, CV_8UC4, pixels);
        cv::Mat flippedImage;
        cv::flip(image, flippedImage, 0);
        cv::Mat rgbImage;
        cv::cvtColor(flippedImage, rgbImage, cv::COLOR_RGBA2BGR);
        delete[] pixels;

        return rgbImage;
    }

    void samplePointsAlongLine(const cv::Point& p1, const cv::Point& p2, int numSamples = 100, 
        int windowWidth = 800, int windowHeight = 600, const cv::Mat& waterDepth = cv::Mat(), int timestamp=0) {
        // scale p1 and p2 to be within the terrain image dimensions
        cv::Point scaledP1(
            static_cast<int>(p1.x * terrainImage.cols / windowWidth),
            static_cast<int>(p1.y * terrainImage.rows / windowHeight)
        );
        cv::Point scaledP2(
            static_cast<int>(p2.x * terrainImage.cols / windowWidth),
            static_cast<int>(p2.y * terrainImage.rows / windowHeight)
        );

        std::vector<float> terrainSampledValues(numSamples);
        std::vector<float> waterSampledValues(numSamples);

        for (int i = 0; i < numSamples; ++i) {
            float alpha = static_cast<float>(i) / (numSamples - 1);
            int x = static_cast<int>(scaledP1.x * (1 - alpha) + scaledP2.x * alpha);
            int y = static_cast<int>(scaledP1.y * (1 - alpha) + scaledP2.y * alpha);
            int wx = static_cast<int>(p1.x * (1 - alpha) + p2.x * alpha);
            int wy = static_cast<int>(p1.y * (1 - alpha) + p2.y * alpha);

            terrainSampledValues[i] = (x >= 0 && x < terrainImage.cols && y >= 0 && y < terrainImage.rows)
                   ? terrainImage.at<uint8_t>(y, x) / 255.0f
                   : 0;
            waterSampledValues[i] = (wx >= 0 && wx < waterDepth.cols && wy >= 0 && wy < waterDepth.rows)
                   ? waterDepth.at<float>(wy, wx) * 2.55f
                   : 0.0f;

            // if (waterSampledValues[i] < terrainSampledValues[i]) {
            //     waterSampledValues[i] = terrainSampledValues[i];
            // }
        }

        std::cout << "Terrain sampled values: ";
        for (const auto& val : terrainSampledValues) std::cout << val << " ";
        std::cout << std::endl;

        std::cout << "Water sampled values: ";
        for (const auto& val : waterSampledValues) std::cout << val << " ";
        std::cout << std::endl;

        std::ofstream outFile("sampled_values.txt", std::ios::app);
        if (outFile.is_open()) {
            if (timestamp == 0) {
            outFile << "Terrain sampled values:\n";
            for (const auto& val : terrainSampledValues) {
                outFile << val << " ";
            }
            outFile << "\n";
            }

            outFile << "Water sampled values (timestamp " << timestamp << "):\n";
            for (const auto& val : waterSampledValues) {
            outFile << val << " ";
            }
            outFile << "\n";

            outFile.close();
            std::cout << "Sampled values appended to sampled_values.txt" << std::endl;
        } else {
            std::cerr << "Failed to open file for appending sampled values." << std::endl;
        }
    }
};

#endif // VISUALISATION_HPP