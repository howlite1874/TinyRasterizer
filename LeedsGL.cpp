#include "LeedsGL.h"
#include <vector>
#include <math.h>
#include <qmath.h>


using namespace std;
LeedsGL::LeedsGL()
{
    //TODO: Initialize all variables that should be initialized!
    inputQueue.clear();
    transformedQueue.clear();
    primitivesQueue.clear();
    clippedPrimitivesQueue.clear();
    fragmentQueue.clear();
    rasterizedPointSize = 0;
    rasterizedLineWidth = 0;

}

LeedsGL::~LeedsGL()
{
    //TODO: Free any resources you allocate.

}

void LeedsGL::clearFramebufferColor()
{
    for (size_t i=0;i<frameBuffer.height;i++)
    {
        for(size_t j=0;j<frameBuffer.width;j++)
        {
            frameBuffer.block[i*frameBuffer.width+j]=bufferClearColor;
        }
    }
}
void LeedsGL::clearDepth(){
    for(size_t i = 0;i < depthBuffer.height;i++)
    {
        for(size_t j = 0;j < depthBuffer.width;j++)
        {
            auto & rgba = depthBuffer[static_cast<int>(i)][static_cast<int>(j)];
            auto depth = reinterpret_cast<float*>(&rgba);
            *depth = 1.f;
        }
    }
}


void LeedsGL::clear(byte mask)
{
    //TODO: Clear the buffer requested using the provided mask.
    //Look at the usage of this function on LeedsGLRenderWidget to
    //know what to expect of the "mask" parameter. Each bit should
    //tell which buffer to clear.


    // 00000001 CLOLOR MASK -> 1
    // 00000010 DEPTH MASK -> 2
    // 00000011 color | depth -> 3

    if( (mask & COLORMASK) == COLORMASK )
    {
        clearFramebufferColor();
    }
    if( (mask & DEPTHMASK) == DEPTHMASK)
    {
        clearDepth();
    }
}

void LeedsGL::setUniform(const string& name,const bool value)
{
    //TODO: set an uniform bool with given name value.
    //uniform variables should decide how your shading happens
    if(name=="texturingEnabled"){
        texturingEnabled=value;
    }

    if(name=="textureModulationEnabled"){
        textureModulationEnabled=value;
    }

    if(name=="UVColorDebug"){
        UVColourDebug=value;
    }

    if(name=="lightingEnabled"){
        lightingEnabled=value;
    }


}

void LeedsGL::setUniform(const string &name, const Matrix4 &mat)
{
    //TODO: set an uniform matrix with given name.
    //uniform variables should decide how your shading happens

    if (name == "viewportMatrix")
    {
        viewPortMatrix = mat;
    }
    else if (name == "projectionMatrix")
    {
        projectionMatrix = mat;
    }
    else if (name == "modelviewMatrix")
    {
        modelviewMatrix = mat;
    }
    else if (name=="lightMatrix")
    {
        lightMatrix = mat;
    }

}

void LeedsGL::setUniform(const string &name, const RGBAValueF &col)
{
    //TODO: set an uniform colour with given name.
    //uniform variables should decide how your shading happens
    if(name=="emissiveMaterial"){
        emissiveMaterial=col;
    }

    if(name=="ambientMaterial"){
        ambientMaterial=col;
    }

    if(name=="diffuseMaterial"){
        diffuseMaterial=col;
    }

    if(name=="specularMaterial"){
        specularMaterial=col;
    }

    if(name=="lightColour"){
        lightColour=col;
    }

}

void LeedsGL::setUniform(const std::string &name, const float val)
{
    //TODO: set an uniform float with given name.
    //uniform variables should decide how your shading happens
    if(name=="shininessMaterial"){
        shininessMaterial=val;
    }

}

void LeedsGL::setUniform(const std::string& name, const Homogeneous4& pos)
{
    //TODO: set an uniform position with given name.
    //uniform variables should decide how your shading happens
    if(name=="lightPosition"){
        lightPosition=pos;
    }

}

void LeedsGL::clearColor(const RGBAValueF &col)
{

    //TODO: set a specific pipeline value, the color used to clear the buffer.
    bufferClearColor=col;

}

void LeedsGL::resizeBuffers(unsigned const int width, unsigned const int height)
{
    //TODO: Implement how you resize your depth buffer
    frameBuffer.Resize(width, height);
    depthBuffer.Resize(width, height);
}

Matrix4 LeedsGLUtils::calculateViewportMatrix(float cx, float cy, float width, float height)
{
    //TODO: return a matrix that given parameters cx,cy (center of your viewport) width height
    //performs the viewport transformation. NDCS -> DCS.
    Matrix4 viewPort;
    viewPort.SetIdentity();
    viewPort.coordinates[0][0]=-width/2;
    viewPort.coordinates[1][1]=-height/2;
    viewPort.coordinates[0][3]=cx+width/2;
    viewPort.coordinates[1][3]=cy+height/2;
    return viewPort;
}

Matrix4 LeedsGLUtils::calculateProjectionOrtho(float left, float right, float bottom, float top, float near, float far)
{
    //TODO: return a Ortographic projection matrix, with the parameters above.
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    Matrix4 ortho;
    ortho.SetIdentity();

    ortho[0][0] = -2.f/(right - left);
    ortho[1][1] = 2.f/(bottom - top) ;
    ortho[2][2] = 1.f/(far - near);

    ortho[0][3] = (left + right) / (left - right);
    ortho[1][3] = (bottom + top) / (bottom - top);
    ortho[2][3] = (near+far)/(near - far);

    return ortho;
}

Matrix4 LeedsGLUtils::calculateProjectionFrustum(float left, float right, float bottom, float top, float near, float far)
{
    //TODO: return a Perspective projection matrix, with the parameters above.
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    Matrix4 Frustum;

    LEFT=-743;
    RIGHT=743;
    BOTTOM=-600;
    TOP=600;

    Frustum[0][0]=2.f * near / (right - left);
    Frustum[1][1]=2.f* near /(top-bottom) ;

    Frustum[0][2] = (left + right) / (left-right);
    Frustum[1][2] = (bottom + top) / (bottom-top);
    Frustum[2][2] = (far + near)/(near-far);

    Frustum[2][3] = -(2 *far * near) / (near - far);

    Frustum[3][2] = 1;


    return Frustum;
}


void LeedsGL::texImage2D(RGBAImage const *textureImage)
{
    //TODO: set in your pipeline which texture should be used to render.
    //Parameter is a pointer to the texture, be aware of how it is stored, and ownership of the resources.
    enabledTexture=textureImage;
}

void LeedsGL::enable(const std::byte function)
{
    //TODO: enables a pipeline function described by the byte parameter.
    if(function==DEPTHTEST)
    {
        //initial value is false
        depthTestEnabled=true;
        if(depthBuffer.width!=frameBuffer.width&&depthBuffer.height!=frameBuffer.height)
        {
            depthBuffer.Resize(frameBuffer.width,frameBuffer.height);
        }
    }
    if(function==PERSPECTIVE)
    {
        perspective=true;
    }

}

void LeedsGL::disable(const std::byte function)
{
    //TODO: disables a pipeline function described by the byte parameter.
    if(function==std::byte{3}){
        depthTestEnabled=false;
    }
    if(function==std::byte{4}){
        perspective=false;
    }

}

void LeedsGL::lineWidth(const float width)
{
    //TODO: Set a variable that describes what is the width in pixels of the line to be rasterized.
    rasterizedLineWidth=width;
}

void LeedsGL::pointSize(const float size)
{
    //TODO: Set a variable that describes what is the size in pixels of the point to be rasterized.
    rasterizedPointSize=size;
}

void LeedsGL::drawArrays(const std::vector<Homogeneous4>& vertices, const std::vector<Homogeneous4>& normals, const std::vector<Cartesian3>& textureCoordinates, const std::vector<RGBAValueF>& colors,std::byte mode)
{
    //Calls the whole pipeline, step by step.

    inputAssembly(vertices,normals,textureCoordinates,colors,inputQueue);
    transformVertices(inputQueue,transformedQueue);
    primitiveAssembly(transformedQueue,mode,primitivesQueue);
    clipAndCull(primitivesQueue,mode,clippedPrimitivesQueue);
    rasterisePrimitives(clippedPrimitivesQueue,mode,fragmentQueue);
    processFragments(fragmentQueue);

    inputQueue.clear();
    transformedQueue.clear();
    primitivesQueue.clear();
    clippedPrimitivesQueue.clear();
    fragmentQueue.clear();


}

void LeedsGL::inputAssembly(const std::vector<Homogeneous4> &vertices, const std::vector<Homogeneous4> &normals, const std::vector<Cartesian3> &textureCoordinates, const std::vector<RGBAValueF> &colors, std::vector<InputVertex> &result)
{
    //    TODO: Check how the input is passed to drawArrays.
    //    This function should combine this disjoint information into a series of InputVertex to be processed
    //    by the next step.
    if(vertices.size()>0){
        for (unsigned int i=0;i<vertices.size();i++) {
            InputVertex tmp(vertices[i]);
            result.emplace_back(tmp);
        }
        for(unsigned int i=0;i<normals.size();i++){
            result[i].normal=normals[i];
        }
        for(unsigned int i=0;i<textureCoordinates.size();i++){
            result[i].texcoord=textureCoordinates[i];
        }
        for(unsigned int i=0;i<colors.size();i++){
            result[i].color=colors[i];
        }

    }

}

void LeedsGL::transformVertices(std::vector<InputVertex> &vertices, std::vector<TransformedVertex>& result)
{
    //TODO: Transform the input vertices using the matrices set in LeedsGLRenderWidget.
    //Also pass all the necessary information to the next steps of the pipeline.
    //You should check the slides to decide which is the appropriate coordinate system to transform them to.
    for (auto& vert : vertices)
    {
        TransformedVertex transVert;
        transVert.position = projectionMatrix * modelviewMatrix * vert.position;
        transVert.color=vert.color;
        transVert.normal=vert.normal.Vector();

        transVert.texcoord=vert.texcoord;
        result.emplace_back(transVert);

    }

}

void LeedsGL::primitiveAssembly(std::vector<TransformedVertex> &vertices,std::byte mode, std::vector<Primitive>& result)
{
    //TODO: Assemble the vertices into a primitive according to the selected mode.
    unsigned int count;
    if(mode==POINTS){
        count=0;
        while(count<vertices.size()){
            for(auto& point:vertices){
                Primitive pri;
                pri.type=mode;
                pri.transformedVertices.emplace_back(point);

                result.emplace_back(pri);
                count++;

            }
        }
    }
    if(mode==LINES){
        count=0;
        while(count<vertices.size()/2){
            Primitive pri;
            pri.type=mode;
            pri.transformedVertices.emplace_back(vertices[count*2]);
            pri.transformedVertices.emplace_back(vertices[count*2+1]);
            result.emplace_back(pri);
            count++;

        }
    }
    if(mode==TRIANGLES){
        count=0;
        while(count<vertices.size()/3){
            Primitive pri;
            pri.type=mode;
            pri.transformedVertices.emplace_back(vertices[count*3]);
            pri.transformedVertices.emplace_back(vertices[count*3+1]);
            pri.transformedVertices.emplace_back(vertices[count*3+2]);
            result.emplace_back(pri);
            count++;
        }
    }
}

//if the point is inside the triangle
OutCode LeedsGL::ComputeOutCode(double x, double y){
    OutCode code;
    code = INSIDECODE;          // initialised as being inside of clip window

    if (x < LeedsGLUtils::LEFT)           // to the left of clip window
        code |= LEFTCODE;
    else if (x >LeedsGLUtils::RIGHT)      // to the right of clip window
        code |= RIGHTCODE;
    if (y < LeedsGLUtils::BOTTOM)           // below the clip window
        code |= BOTTOMCODE;
    else if (y > LeedsGLUtils::TOP)      // above the clip window
        code |= TOPCODE;

    return code;
}

// Get intersection position of two lines.
Line LeedsGL::GetCrossVertex(TransformedVertex v0, TransformedVertex v1){
    Line res;
    TransformedVertex vout, vin;

    OutCode outcode0 = ComputeOutCode(v0.position.x, v0.position.y);
    OutCode outcode1 = ComputeOutCode(v1.position.x, v1.position.y);

    while (true){
        if (!(outcode0 | outcode1)){//if true two points are inside the plane,so exit the circulation
            return Line(v0,v1);
        }
        else if (outcode0 & outcode1){// if true two points are all outside the plane,so exit the circulation
            return res;
        }
        else{
            //find point which is outside the plane
            OutCode outcodeOut = outcode0 ? outcode0 : outcode1;//outcodeOut is the code of point which is outside the plane
            vin = (outcode0 == outcodeOut) ? v1 : v0;
            vout = (outcode0 == outcodeOut) ? v0 : v1;
            float x1 = vin.position.x,
                    y1 = vin.position.y,
                    z1 = vin.position.z,
                    w1 = vin.position.w;
            float x2 = vout.position.x,
                    y2 = vout.position.y,
                    z2 = vout.position.z,
                    w2 = vout.position.w;
            TransformedVertex crossPoint;

            // find the intersection
            if (outcodeOut & TOPCODE){
                if (y2 == y1)
                    return res;
                crossPoint.position.x = x1 + (x2 - x1) * (LeedsGLUtils::TOP - y1) / (y2 - y1);
                crossPoint.position.y = LeedsGLUtils::TOP;
            }
            else if (outcodeOut & BOTTOMCODE){
                if (y2 == y1)
                    return res;
                crossPoint.position.x = x1 + (x2 - x1) * (LeedsGLUtils::BOTTOM - y1) / (y2 - y1);
                crossPoint.position.y = LeedsGLUtils::BOTTOM;
            }
            else if (outcodeOut & LEFTCODE){
                if (x2 == x1)
                    return res;
                crossPoint.position.y = y1 + (y2 - y1) * (LeedsGLUtils::RIGHT - x1) / (x2 - x1);
                crossPoint.position.x = LeedsGLUtils::RIGHT;
            }
            else if (outcodeOut & RIGHTCODE){
                if (x2 == x1)
                    return res;
                crossPoint.position.y = y1 + (y2 - y1) * (LeedsGLUtils::LEFT - x1) / (x2 - x1);
                crossPoint.position.x = LeedsGLUtils::LEFT;
            }

            float dy = crossPoint.position.y - v0.position.y;
            float t = dy / (v1.position.y - v0.position.y);
            crossPoint.texcoord = (1-t)*v0.texcoord+ t*v1.texcoord;
            crossPoint.normal = (1-t)*v0.normal+t * v1.normal;
            crossPoint.color = (1-t)*v0.color+ t * v1.color;


            //in case the other point is outside the plane
            if (outcodeOut == outcode0) {
                v0 = crossPoint;
                outcode0 = ComputeOutCode(v0.position.x, v0.position.y);
            }
            else {
                v1 = crossPoint;
                outcode1 = ComputeOutCode(v1.position.x, v1.position.y);
            }
        }
    }
}


void LeedsGL::clipAndCull(std::vector<Primitive>& primitives,std::byte mode, std::vector<Primitive>& result)
{
    //TODO: Implement clipping and culling. Should have a different behavior for each type of primitive.
    //Pay attention to what type of projection you are using, as your clipping planes will be different.
    //If you choose to skip this step as it is one of your last tasks, just return all the same primitives.
    float minX = frameBuffer.width;
    float maxX = 0.0;
    float minY = frameBuffer.height;
    float maxY = 0.0;

    if(mode==POINTS){
        for(auto& p:primitives){
            if (p.transformedVertices[0].position.x < minX) minX = p.transformedVertices[0].position.x;
            if (p.transformedVertices[0].position.x > maxX) maxX = p.transformedVertices[0].position.x;
            if (p.transformedVertices[0].position.y < minY) minY = p.transformedVertices[0].position.y;
            if (p.transformedVertices[0].position.y > maxY) maxY = p.transformedVertices[0].position.y;
            result.emplace_back(p);
        }
    }
    else if(mode==LINES){
        for(auto& p:primitives){
            if (p.transformedVertices[0].position.x < minX) minX = p.transformedVertices[0].position.x;
            if (p.transformedVertices[1].position.x < minX) minX = p.transformedVertices[1].position.x;

            if (p.transformedVertices[0].position.x > maxX) maxX = p.transformedVertices[0].position.x;
            if (p.transformedVertices[1].position.x > maxX) maxX = p.transformedVertices[1].position.x;

            if (p.transformedVertices[0].position.y < minY) minY = p.transformedVertices[0].position.y;
            if (p.transformedVertices[1].position.y < minY) minY = p.transformedVertices[1].position.y;

            if (p.transformedVertices[0].position.y > maxY) maxY = p.transformedVertices[0].position.y;
            if (p.transformedVertices[1].position.y > maxY) maxY = p.transformedVertices[1].position.y;
            result.emplace_back(p);
        }
    }
    else if(mode==TRIANGLES){
        Boundary leftBoundary(0, 1, -frameBuffer.width);                     // x = 0
        Boundary rightBoundary(0, 1, frameBuffer.width);   // x = frameBuffer.width
        Boundary bottomBoundary(0, 1, -frameBuffer.height);                   // y = 0
        Boundary topBoundary(0, 1, frameBuffer.height);    // y = frameBuffer.height
        for(auto& p:primitives){
            //cull
            // Calculate face normal of the triangle
            Cartesian3 edge1 = (p.transformedVertices[1].position - p.transformedVertices[0].position).Vector();
            Cartesian3 edge2 = (p.transformedVertices[2].position - p.transformedVertices[0].position).Vector();
            Cartesian3 faceNormal = edge1.cross(edge2);

            Cartesian3 viewDirection(0, 0, 1);

            if(faceNormal.dot(viewDirection) <= 0) {
                continue;  // Skip this triangle
            }

            //clip
            std::vector<TransformedVertex> vertices = {
                        p.transformedVertices[0],
                        p.transformedVertices[1],
                        p.transformedVertices[2]
                    };

            for (auto& boundary : {leftBoundary, rightBoundary, bottomBoundary, topBoundary}) {
                vertices = SutherlandHodgmanClip(vertices, boundary);
            }

            // Depending on the number of vertices returned by clipping,
            // you might end up with triangles, quadrilaterals or more complex polygons.
            if(vertices.size() == 3) {
                // It remains a triangle. You can push it to your result or process it directly.
                Primitive newTri;
                newTri.transformedVertices.push_back(vertices[0]);
                newTri.transformedVertices.push_back(vertices[1]);
                newTri.transformedVertices.push_back(vertices[2]);
                result.emplace_back(newTri);
            } else if(vertices.size() == 4) {
                // It became a quadrilateral. Split it into two triangles.
                Primitive tri1, tri2;
                tri1.transformedVertices.push_back(vertices[0]);
                tri1.transformedVertices.push_back(vertices[1]);
                tri1.transformedVertices.push_back(vertices[2]);


                tri2.transformedVertices.push_back(vertices[0]);
                tri2.transformedVertices.push_back(vertices[2]);
                tri2.transformedVertices.push_back(vertices[3]);

                result.emplace_back(tri1);
                result.emplace_back(tri2);
            }


        }
    }
}

bool LeedsGL::isInside(const TransformedVertex& vertex, const Boundary& boundary) {
    return boundary.A * vertex.position.x + boundary.B * vertex.position.y + boundary.C >= 0;
}

std::vector<TransformedVertex> LeedsGL::SutherlandHodgmanClip(const std::vector<TransformedVertex>& inputVertices, const Boundary& boundary) {
    std::vector<TransformedVertex> outputVertices;

    for(size_t i = 0; i < inputVertices.size(); ++i) {
        TransformedVertex currentVertex = inputVertices[i];
        TransformedVertex nextVertex = inputVertices[(i + 1) % inputVertices.size()];

        // Check if currentVertex and nextVertex are inside the boundary
        if(isInside(currentVertex, boundary)) {
            if(isInside(nextVertex, boundary)) {
                outputVertices.push_back(nextVertex);
            } else {
                TransformedVertex intersection = computeIntersection(currentVertex, nextVertex, boundary);
                outputVertices.push_back(intersection);
            }
        } else if(isInside(nextVertex, boundary)) {
            TransformedVertex intersection = computeIntersection(currentVertex, nextVertex, boundary);
            outputVertices.push_back(intersection);
            outputVertices.push_back(nextVertex);
        }
    }

    return outputVertices;
}

TransformedVertex LeedsGL::computeIntersection(const TransformedVertex& start, const TransformedVertex& end, const Boundary& boundary) {
    // Direction vector
    Cartesian3 d = (end.position - start.position).Vector();

    // Compute t
    float denominator = boundary.A * d.x + boundary.B * d.y;

    // Ensure we're not dividing by zero (i.e., line segment is not parallel to boundary)
    if(fabs(denominator) < 1e-6) { // some small threshold
        // Handle this case - maybe return start or end or some sentinel value
    }

    float t = -(boundary.A * start.position.x + boundary.B * start.position.y + boundary.C) / denominator;

    // Compute intersection point
    TransformedVertex intersectionVertex;
    intersectionVertex.position = start.position + t * d;

    // If you have other attributes (color, texture coordinates, etc.), interpolate them here based on t.
    // ...

    return intersectionVertex;
}



void LeedsGL::rasterisePrimitives(std::vector<Primitive> &primitives, std::byte mode, std::vector <Fragment>& results)
{
    //TODO: Generate a list of fragments according to what mode is chosen. Should call the "rasterise X" functions
    if(mode==POINTS){
        for (auto& pri:primitives) {
            pri.transformedVertices[0].position=viewPortMatrix * pri.transformedVertices[0].position.Point();
            //world -> NDC[0,1] -> x,y

            float zFar = 10.f;
            float zNear = 0.1f;

            //save non-linear-z
            pri.transformedVertices[0].position.z =  ( (zFar - zNear) *  pri.transformedVertices[0].position .z + (zFar+zNear) ) * 0.5f;
            rasterisePoint(pri,results);
        }
    }
    if(mode==LINES){
        for (auto& pri:primitives) {
            pri.transformedVertices[0].position=viewPortMatrix * pri.transformedVertices[0].position.Point();
            pri.transformedVertices[1].position=viewPortMatrix * pri.transformedVertices[1].position.Point();

            float zFar = 10.f;
            float zNear = 0.1f;

            //save non-linear-z
            pri.transformedVertices[0].position.z =  ( (zFar - zNear) *  pri.transformedVertices[0].position .z + (zFar+zNear) ) * 0.5f;
            pri.transformedVertices[1].position.z =  ( (zFar - zNear) *  pri.transformedVertices[1].position .z + (zFar+zNear) ) * 0.5f;

            rasteriseLine(pri,results);

        }
    }
    if(mode==TRIANGLES){
        for(auto& pri:primitives){
            pri.transformedVertices[0].position=viewPortMatrix * pri.transformedVertices[0].position.Point();
            pri.transformedVertices[1].position=viewPortMatrix * pri.transformedVertices[1].position.Point();
            pri.transformedVertices[2].position=viewPortMatrix * pri.transformedVertices[2].position.Point();

            float zFar = 10.f;
            float zNear = 0.1f;

            //save non-linear-z
            pri.transformedVertices[0].position.z =  ( (zFar - zNear) *  pri.transformedVertices[0].position .z + (zFar+zNear) ) * 0.5f;
            pri.transformedVertices[1].position.z =  ( (zFar - zNear) *  pri.transformedVertices[1].position .z + (zFar+zNear) ) * 0.5f;
            pri.transformedVertices[2].position.z =  ( (zFar - zNear) *  pri.transformedVertices[2].position .z + (zFar+zNear) ) * 0.5f;

            rasteriseTriangle(pri,results);
        }
    }
}

bool LeedsGL::depthTest(float x,float y,float z){
    //single channel did have enough precision.
    //using 4 channels to represtent single depth(float)
    //single float is 32bit, rgab is 32bit.

    auto & rgba = depthBuffer[static_cast<int>(x)][static_cast<int>(y)];
    float& depth = reinterpret_cast<float&>(rgba);

    if(depth >z) return true;
    return false;
}

void LeedsGL::rasterisePoint(const Primitive &point,std::vector<Fragment>& output)
{
    //TODO: Rasterise a point, according to the pointSize.

    //initialize one point
    Fragment tmp;

    tmp.normal=point.transformedVertices[0].normal;
    tmp.position=point.transformedVertices[0].position.Point();
    tmp.texCoord=point.transformedVertices[0].texcoord;
    tmp.color=point.transformedVertices[0].color;

    int startX = point.transformedVertices[0].position.x -  rasterizedPointSize / 2;
    int startY = point.transformedVertices[0].position.y -  rasterizedPointSize / 2;

    if(rasterizedPointSize>0){
        for (unsigned int i = 0; i < rasterizedPointSize; i++) {
            for (unsigned int j = 0; j < rasterizedPointSize; j++) {
                tmp.row = startY + j;
                tmp.col = startX + i;
                if(depthTestEnabled)
                {
                    if(depthTest(tmp.row,tmp.col,tmp.position.z)){

                        auto & rgba = depthBuffer[tmp.row][tmp.col];
                        auto depth = reinterpret_cast<float*>(&rgba);
                        *depth = tmp.position.z;

                        fragmentQueue.push_back(tmp);
                    }
                }
                else
                {
                    fragmentQueue.push_back(tmp);
                }
            }
        }

    }
    else{
        tmp.row= point.transformedVertices[0].position.y;
        tmp.col= point.transformedVertices[0].position.x;
        if(depthTestEnabled)
        {
            if(depthTest(tmp.row,tmp.col,tmp.position.z)){

                auto & rgba = depthBuffer[tmp.row][tmp.col];
                auto depth = reinterpret_cast<float*>(&rgba);
                *depth = tmp.position.z;

                fragmentQueue.push_back(tmp);
            }
        }
        else
        {
            fragmentQueue.push_back(tmp);
        }
    }


}

void LeedsGL::rasteriseLine(const Primitive &line, std::vector<Fragment> &output)
{
    TransformedVertex vertBegin = line.transformedVertices[0];
    TransformedVertex vertEnd = line.transformedVertices[1];

    Cartesian3 dpos = vertEnd.position.Point() - vertBegin.position.Point();
    float length = sqrt(dpos.x * dpos.x + dpos.y * dpos.y);
    float dx = dpos.x / length;
    float dy = dpos.y / length;

    for (float t = 0.f; t <= 1.0f; t += 0.001f) {
        TransformedVertex tmp;
        tmp.position = vertBegin.position + dpos * t;
        tmp.color.red = (1-t) * vertBegin.color.red + t * vertEnd.color.red;
        tmp.color.green = (1-t) * vertBegin.color.green + t * vertEnd.color.green;
        tmp.color.blue = (1-t) * vertBegin.color.blue + t * vertEnd.color.blue;
        tmp.color.alpha = (1-t) * vertBegin.color.alpha + t * vertEnd.color.alpha;

        for (unsigned int i = 0; i < rasterizedLineWidth; i++) {
            Fragment frag;
            float offset = (float)i - rasterizedLineWidth / 2.0f;

            // Offset based on the slope of the line
            frag.row = static_cast<int>((tmp.position.y + offset * dx));
            frag.col = static_cast<int>((tmp.position.x - offset * dy));
            frag.color = tmp.color;

            if (depthTestEnabled) {
                if (depthTest(frag.col, frag.row, tmp.position.z)) {
                    auto & rgba = depthBuffer[frag.row][frag.col];
                    auto depth = reinterpret_cast<float*>(&rgba);
                    *depth = tmp.position.z;
                    fragmentQueue.push_back(frag);
                }
            } else {
                fragmentQueue.push_back(frag);
            }
        }
    }
}


float LeedsGLUtils::distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p)
{
    //assumes n is normalized
    return n.dot(r) - n.dot(p);
}
void LeedsGL::rasteriseTriangle(const Primitive &triangle, std::vector<Fragment> &output)
{
    TransformedVertex vertex0 = triangle.transformedVertices[0];
    TransformedVertex vertex1 = triangle.transformedVertices[1];
    TransformedVertex vertex2 = triangle.transformedVertices[2];

    // compute a bounding box that starts inverted to frame size
    // clipping will happen in the raster loop proper
    float minX = frameBuffer.width, maxX = 0.0;
    float minY = frameBuffer.height, maxY = 0.0;

    // test against all vertices
    if (vertex0.position.x < minX) minX = vertex0.position.x;
    if (vertex0.position.x > maxX) maxX = vertex0.position.x;
    if (vertex0.position.y < minY) minY = vertex0.position.y;
    if (vertex0.position.y > maxY) maxY = vertex0.position.y;

    if (vertex1.position.x < minX) minX = vertex1.position.x;
    if (vertex1.position.x > maxX) maxX = vertex1.position.x;
    if (vertex1.position.y < minY) minY = vertex1.position.y;
    if (vertex1.position.y > maxY) maxY = vertex1.position.y;

    if (vertex2.position.x < minX) minX = vertex2.position.x;
    if (vertex2.position.x > maxX) maxX = vertex2.position.x;
    if (vertex2.position.y < minY) minY = vertex2.position.y;
    if (vertex2.position.y > maxY) maxY = vertex2.position.y;

    Cartesian3 v0 = Cartesian3(vertex0.position.x,vertex0.position.y,0);
    Cartesian3 v1 = Cartesian3(vertex1.position.x,vertex1.position.y,0);
    Cartesian3 v2 = Cartesian3(vertex2.position.x,vertex2.position.y,0);

    Cartesian3 v0v1 = v1-v0;
    Cartesian3 n_v0v1 = Cartesian3(-v0v1.y,v0v1.x,0);
    Cartesian3 v1v2 = v2-v1;
    Cartesian3 n_v1v2 = Cartesian3(-v1v2.y,v1v2.x,0);
    Cartesian3 v2v0 = v0-v2;
    Cartesian3 n_v2v0 = Cartesian3(-v2v0.y,v2v0.x,0);

    float dAlpha =  LeedsGLUtils::distancePointLine(v0, n_v1v2,v1);
    float dBeta = LeedsGLUtils::distancePointLine(v1, n_v2v0,v2);
    float dGamma = LeedsGLUtils::distancePointLine(v2,n_v0v1,v0);

    if (abs(dAlpha-0)<std::numeric_limits<float>::epsilon() ||
            abs(dBeta-0)<std::numeric_limits<float>::epsilon() ||
            abs(dGamma-0)<std::numeric_limits<float>::epsilon())
        return;

    // create a fragment for reuse
    Fragment rasterFragment;

    // loop through the pixels in the bounding box
    for (rasterFragment.row = int(minY); rasterFragment.row <= maxY; rasterFragment.row++)
    { // per row
        // this is here so that clipping works correctly
        if (rasterFragment.row < 0) continue;
        if (rasterFragment.row >= int(frameBuffer.height)) continue;
        for (rasterFragment.col = int(minX); rasterFragment.col <= maxX; rasterFragment.col++)
        { // per pixel
            // this is also for correct clipping
            if (rasterFragment.col < 0) continue;
            if (rasterFragment.col >= int(frameBuffer.width)) continue;

            // the pixel in cartesian format
            Cartesian3 pixel(rasterFragment.col+0.5f, rasterFragment.row+0.5f, 0.0f);

            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = LeedsGLUtils::distancePointLine(pixel,n_v1v2,v1) / dAlpha;
            float beta = LeedsGLUtils::distancePointLine(pixel,n_v2v0,v2)/dBeta;
            float gamma = LeedsGLUtils::distancePointLine(pixel,n_v0v1,v0)/dGamma;

            // now perform the half-plane test
            if ((alpha < 0.0f) || (beta < 0.0f) || (gamma < 0.0f))
                continue;

            //TODO: use the computed baricentric coordinates to interpolate all of the necessary properties for rendering.
            //Be aware of making sure they are perspective correct.

            rasterFragment.color = alpha * vertex0.color + beta * vertex1.color + gamma * vertex2.color;
            rasterFragment.texCoord =  (alpha * vertex0.texcoord + beta * vertex1.texcoord + gamma * vertex2.texcoord);
            rasterFragment.normal =alpha * vertex0.normal + beta * vertex1.normal + gamma * vertex2.normal;

            auto vertex = alpha * vertex0.position + beta * vertex1.position + gamma * vertex2.position;
            if(depthTestEnabled)
            {
                if(depthTest(rasterFragment.row,rasterFragment.col,vertex.z)){

                    auto & rgba = depthBuffer[rasterFragment.row][rasterFragment.col];
                    auto depth = reinterpret_cast<float*>(&rgba);
                    *depth = vertex.z;

                    fragmentQueue.push_back(rasterFragment);
                }
            }
            else
            {
                fragmentQueue.push_back(rasterFragment);
            }
        } // per pixel
    } // per row
}

void LeedsGL::processFragments(std::vector<Fragment> &fragments)
{
    //TODO: Process all of the fragments, shading according to the uniform properties.
    //Depth test should go here. We don't explicitly have a pre or post fragment stage.
    //Consider the "shading" as the fragment stage. Decide if the depth test should go before or after, and justify.
#pragma omp parallel
    for(auto& frag: fragments){
        if(frag.row<=frameBuffer.height&&frag.col<=frameBuffer.width && frag.row >= 0 && frag.col >= 0){
            RGBAValueF light = RGBAValueF(1.f, 1.f, 1.f, 1.f);
            if(lightingEnabled) light = CalculateLighting(frag.normal,frag.position,this->emissiveMaterial, this->ambientMaterial,this->diffuseMaterial, this->specularMaterial, this->shininessMaterial);
            if(this->texturingEnabled)
            {
                float u = min(int(frag.texCoord.x * (enabledTexture->width )), int(enabledTexture->width));
                float v = min(int(frag.texCoord.y * (enabledTexture->height )), int(enabledTexture->height));

                //bilinear interpolation
                int s = static_cast<int>(u);
                int t = static_cast<int>(v);
                float sParm = u - s;
                float tParm = v - t;

                //grab nearest texture colors
                RGBAValue color00 = (*enabledTexture)[s % enabledTexture->width][t % enabledTexture->height];
                RGBAValue color01 = (*enabledTexture)[s % enabledTexture->width][(t + 1) % enabledTexture->height];
                RGBAValue color10 = (*enabledTexture)[(s + 1) % enabledTexture->width][t % enabledTexture->height];
                RGBAValue color11 = (*enabledTexture)[(s + 1) % enabledTexture->width][(t + 1) % enabledTexture->height];

                //compute color on edges
                RGBAValue color0 =(1 - sParm) * color00  + sParm * color10;
                RGBAValue color1 = (1 - sParm) *color01  + sParm * color11;

                RGBAValue color = (1 - tParm) * color0  + tParm * color1;

                if(textureModulationEnabled)
                    frameBuffer[frag.row][frag.col] = RGBAValueF(float(color.red)/255.f,float(color.green)/255.f,float(color.blue)/255.f,float(color.alpha)/255.f) * light;
                else
                    frameBuffer[frag.row][frag.col] = RGBAValueF(float(color.red)/255.f,float(color.green)/255.f,float(color.blue)/255.f,float(color.alpha)/255.f) ;
            }
            else {
                frameBuffer[frag.row][frag.col] = frag.color* light;
            }
        }
    }
}

RGBAValueF LeedsGL::CalculateLighting(const Homogeneous4& n_vcs, const Homogeneous4& v_vcs,const RGBAValueF& em, const RGBAValueF& am, const RGBAValueF& diff, const RGBAValueF& spec, float shin)
{    
    if(n_vcs.x == 0.0f && n_vcs.y == 0.0f && n_vcs.z == 0.0f) // we shouldn't try shading if there are no normals
        return RGBAValueF();

    Cartesian3 lightVector;
    Cartesian3 unitNormal = n_vcs.Vector().unit();

    Homogeneous4 lp = lightMatrix * lightPosition;

    if(abs(lp.w - 0) < std::numeric_limits<float>::epsilon())
        lightVector = lp.Vector().unit();
    else //point light
        lightVector = (lp - v_vcs).Vector().unit();
    Cartesian3 eyeVector = perspective? -1*v_vcs.Point(): Cartesian3(0,0,1);
    Cartesian3 bisector = (lightVector + eyeVector).unit();

    RGBAValueF emissive = em;
    RGBAValueF ambient =  am.modulate(ambientMaterial);

    float dDot = unitNormal.dot(lightVector);
    dDot = dDot <0? 0: dDot;

    RGBAValueF diffuse = dDot * diff.modulate(lightColour);

    float sDot = unitNormal.dot(bisector);
    sDot = sDot <0? 0:sDot;
    sDot = pow(sDot,shin);
    //sDot = ((f.shininess+2)/8.0f)*sDot*dDot;
    sDot = dDot>0? sDot : 0;
    sDot = sDot* dDot *(shin+2)/2*float(M_PI);

    Cartesian3 fs = Cartesian3(spec.red,spec.green,spec.blue);
    Cartesian3 air = Cartesian3(1,1,1);
    Cartesian3 a = (air - fs);
    Cartesian3 b = (air + fs);
    Cartesian3 r0 = Cartesian3(a.x/b.x,a.y/b.y,a.z/b.z);
    r0 = Cartesian3(r0.x*r0.x,r0.y*r0.y,r0.z*r0.z);
    Cartesian3 rschlick = r0 + (air-r0) * powf((1-bisector.dot(lightVector)),5);
    RGBAValueF updatedSpecular = RGBAValueF(rschlick.x,rschlick.y,rschlick.z,1);
    RGBAValueF specular =sDot * updatedSpecular.modulate(lightColour);
    return emissive + ambient + diffuse + specular;
}



