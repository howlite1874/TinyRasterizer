#ifndef LEEDSGL_H
#define LEEDSGL_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <cstddef>
#include "Matrix4.h"
#include "RGBAValue.h"
#include "RGBAValueF.h"
#include "RGBAImage.h"

typedef int OutCode;
const int INSIDECODE = 0; // 0000
const int LEFTCODE = 1;   // 0001
const int RIGHTCODE = 2; // 0010
const int BOTTOMCODE = 4; // 0100
const int TOPCODE = 8;    // 1000

namespace LeedsGLUtils
{
    Matrix4 calculateViewportMatrix(float cx, float cy, float width, float height);
    Matrix4 calculateProjectionFrustum(float left, float right,float bottom, float top, float near, float far);
    Matrix4 calculateProjectionOrtho(float left, float right,float bottom, float top, float near, float far);
    float distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p);

    //Frustum
    static float LEFT;
    static float RIGHT;
    static float BOTTOM;
    static float TOP;
}

// class with vertex attributes
struct InputVertex
{
    //TODO: Complete with what information needs to be passed forward.
    Homogeneous4 position;
    RGBAValueF color;
    Cartesian3 texcoord;
    Homogeneous4 normal;

    InputVertex(){
        position = Homogeneous4();
        normal = Homogeneous4();
        texcoord = Cartesian3();
        color = RGBAValueF();
    }

    InputVertex(Homogeneous4 value){
        position = value;
        normal = Homogeneous4();
        texcoord = Cartesian3();
        color = RGBAValueF();

    }

};

struct TransformedVertex
{
    //TODO: Complete with what information needs to be passed forward.
     Homogeneous4 position;
     RGBAValueF color;
     Cartesian3 texcoord;
     Cartesian3 normal;

     TransformedVertex( )
     {
         position = Homogeneous4();
         normal = Cartesian3();
         texcoord = Cartesian3();
         color = RGBAValueF();
     }


};

struct Primitive
{
    //TODO: Complete with what information needs to be passed forward.
    std::byte type;
    std::vector<TransformedVertex> transformedVertices;
    std::vector<unsigned int> elements;

    Primitive()
     {
         type = std::byte{0};
         transformedVertices.clear();
         elements.clear();
     }

};


struct Fragment
{
    int row;
    int col;
    //TODO: Complete with what information needs to be passed forward.
    Cartesian3 position;
    RGBAValueF color;
    Cartesian3 normal;
    Cartesian3 texCoord;


    Fragment(){
        row=0;
        col=0;
        position=Cartesian3();
        color=RGBAValueF();
        normal=Cartesian3();
        texCoord=Cartesian3();
    }

    Fragment(Primitive& point)
    {
        row = 0;
        col = 0;
        if (point.type == std::byte{ 5 })
        {
            position = point.transformedVertices[0].position.Point();
            normal = point.transformedVertices[0].normal;
            texCoord = point.transformedVertices[0].texcoord;
            color = point.transformedVertices[0].color;
        }
    }

};

struct Line
    {
        // membes of a line
        TransformedVertex v1;
        TransformedVertex v2;
        bool isNull = false;

        Line() :isNull(true){}
        Line(const TransformedVertex& v1, const TransformedVertex& v2) :v1(v1), v2(v2) {}
    };

struct Boundary {
    float A, B, C;

    Boundary(float a, float b, float c) : A(a), B(b), C(c) {}
};


class LeedsGL
{
public:
    LeedsGL();
    ~LeedsGL();

    //RENDERING PARAMETERS:
    void setUniform(const std::string& name,const bool value);
    void setUniform(const std::string& name, const Matrix4& mat);
    void setUniform(const std::string& name, const RGBAValueF& col);
    void setUniform(const std::string& name, const Homogeneous4& pos);
    void setUniform(const std::string& name, const float val);

    //PIPELINE CONTROL:
    void clearColor(const RGBAValueF &col);
    void clear(std::byte mask);
    void enable(const std::byte function);
    void disable(const std::byte function);
    void texImage2D(RGBAImage const *textureImage);
    void resizeBuffers(unsigned const int width,unsigned const int height);
    void lineWidth(const float width);
    void pointSize(const float size);

    //MAIN PIPELINE IMPLEMENTATION
    void drawArrays(const std::vector<Homogeneous4>& vertices,
                    const std::vector<Homogeneous4>& normals,
                    const std::vector<Cartesian3>& textureCoordinates,
                    const std::vector<RGBAValueF>& colors,std::byte mode);
    void inputAssembly(const std::vector<Homogeneous4>& vertices,
                       const std::vector<Homogeneous4>& normals,
                       const std::vector<Cartesian3>& textureCoordinates,
                       const std::vector<RGBAValueF>& colors,
                       std::vector<InputVertex>& result);
    void transformVertices(std::vector<InputVertex>& vertices,
                           std::vector<TransformedVertex>& result);
    void primitiveAssembly(std::vector<TransformedVertex>& vertices,
                             std::byte mode,
                           std::vector<Primitive>& result);
    void clipAndCull(std::vector<Primitive>& primitives,
                     std::byte mode,
                     std::vector<Primitive>& result);
    void rasterisePrimitives(std::vector<Primitive>& primitives,
                             std::byte mode,
                             std::vector<Fragment>& result);

    void rasterisePoint(const Primitive& point,
                        std::vector<Fragment>& output);
    void rasteriseLine(const Primitive& line,
                       std::vector<Fragment>& output);
    void rasteriseTriangle(const Primitive& triangle,
                           std::vector<Fragment>& output);
    void processFragments(std::vector<Fragment>& fragments);

    void clearFramebufferColor();

    bool depthTest(float x,float y,float z);

    void clearDepth();

    OutCode ComputeOutCode(double x, double y);

    Line GetCrossVertex(TransformedVertex v0, TransformedVertex v1);

    bool isInside(const TransformedVertex& vertex, const Boundary& boundary) ;

    std::vector<TransformedVertex> SutherlandHodgmanClip(const std::vector<TransformedVertex>& inputVertices, const Boundary& boundary);

    TransformedVertex computeIntersection(const TransformedVertex& start, const TransformedVertex& end, const Boundary& boundary);
    //SHADING.
    RGBAValueF CalculateLighting(const Homogeneous4& n_vcs,
                                 const Homogeneous4& v_vcs,
                                 const RGBAValueF& em,
                                 const RGBAValueF& am,
                                 const RGBAValueF& diff,
                                 const RGBAValueF& spec,
                                 float shin);
    //BUFFERS
    RGBAImage frameBuffer;
    RGBAImage depthBuffer;

    //Masks
    static const std::byte COLORMASK{1};
    static const std::byte DEPTHMASK{2};
    //Function constants
    static const std::byte DEPTHTEST{3};
    static const std::byte PERSPECTIVE{4};
    //Drawing modes

    static const std::byte POINTS{5};
    static const std::byte LINES{6};
    static const std::byte TRIANGLES{7};



private:

//uniform variables
    RGBAImage const* enabledTexture;
    bool texturingEnabled;
    bool textureModulationEnabled;
    bool UVColourDebug;
    bool lightingEnabled;

    Matrix4 lightMatrix;
    Matrix4 viewPortMatrix;
    Matrix4 projectionMatrix;
    Matrix4 modelviewMatrix;

    Homogeneous4 lightPosition;
    RGBAValueF lightColour;

    RGBAValueF emissiveMaterial;
    RGBAValueF ambientMaterial;
    RGBAValueF diffuseMaterial;
    RGBAValueF specularMaterial;
    float shininessMaterial;

//global states
    float rasterizedLineWidth;
    float rasterizedPointSize;
    RGBAValueF bufferClearColor;
    bool depthTestEnabled;
    bool perspective;

//Queues
    std::vector<InputVertex> inputQueue;
    std::vector<TransformedVertex> transformedQueue;
    std::vector<Primitive> primitivesQueue;
    std::vector<Primitive> clippedPrimitivesQueue;
    std::vector<Fragment> fragmentQueue;


};

#endif // LEEDSGL_H
