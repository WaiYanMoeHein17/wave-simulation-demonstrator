#version 330 core

out vec4 FragColor;

in vec2 TexCoord;

// Textures
uniform sampler2D terrain;
uniform sampler2D waterHeightMap;
uniform sampler2D waterJetMap;

// Uniforms
uniform vec3 colorMap[8];
uniform float contourLineFactor;
uniform bool gradientColor;
uniform bool grayscale;

uniform float simulationScale;
uniform float simulationOffset;

uniform bool useWaterTexture;
uniform bool showSimulation;

void main() {
    // https://web.cs.ucdavis.edu/~okreylos/ResDev/SARndbox/
    float corner0 = texture(terrain, TexCoord).r;
    float corner1 = texture(terrain, TexCoord + vec2(1.0 / textureSize(terrain, 0).x, 0.0)).r;
    float corner2 = texture(terrain, TexCoord + vec2(0.0, 1.0 / textureSize(terrain, 0).y)).r;
    float corner3 = texture(terrain, TexCoord + vec2(1.0 / textureSize(terrain, 0).x, 1.0 / textureSize(terrain, 0).y)).r;
    // ^

    float minEl = min(min(corner0, corner1), min(corner2, corner3));
    float maxEl = max(max(corner0, corner1), max(corner2, corner3));

    vec4 baseColor;

    if (floor(maxEl * contourLineFactor) != floor(minEl * contourLineFactor)) {
        baseColor = vec4(0.0, 0.0, 0.0, 1.0);
    } else if (grayscale) {
        baseColor = vec4(corner0, corner0, corner0, 1.0);
    } else if (gradientColor) {
        float gradientIndex = corner0 * 7.0;
        int lowerIndex = int(floor(gradientIndex));
        int upperIndex = int(ceil(gradientIndex));
        float blendFactor = fract(gradientIndex);

        vec3 lowerColor = colorMap[clamp(lowerIndex, 0, 7)];
        vec3 upperColor = colorMap[clamp(upperIndex, 0, 7)];
        vec3 blendedColor = mix(lowerColor, upperColor, blendFactor);

        baseColor = vec4(blendedColor, 1.0);
    } else {
        int index = int(corner0 * 7.0);
        baseColor = vec4(colorMap[clamp(index, 0, 7)], 1.0);
    }

    float waterCanvasValue = texture(waterHeightMap, TexCoord).r * simulationScale + simulationOffset;
    float offset = (corner0 - 0.001) * float(showSimulation);
    // float offset = 0.0;

    if (corner0 < waterCanvasValue + offset) {
        if (useWaterTexture) {
            vec4 waterColor = texture(waterJetMap, TexCoord).rgba;
            float shallowFactor = clamp((waterCanvasValue - corner0) / 0.1, 0.0, 1.0);
            waterColor = mix(baseColor, waterColor, shallowFactor);
            baseColor = mix(baseColor, waterColor, 0.95);
        } else {
            float depthFactor = (waterCanvasValue - corner0) / waterCanvasValue;
            vec4 deepWaterColor = vec4(0.0, 0.0, 0.5, 1.0);
            vec4 shallowWaterColor = vec4(0.0, 1.0, 1.0, 1.0);
            vec4 waterColor = mix(shallowWaterColor, deepWaterColor, depthFactor);

            float shallowFactor = clamp((waterCanvasValue - corner0) / 0.1, 0.0, 1.0);
            waterColor = mix(baseColor, waterColor, shallowFactor);
            baseColor = mix(baseColor, waterColor, 0.95);
        }
    }
    // if(showSimulation) {
    //     float r = texture(waterHeightMap, TexCoord).r;
    //     baseColor = vec4(r,r,r,1.0);
    // }

    FragColor = baseColor;
}
