#ifndef _MY_MESH_
#define _MY_MESH_


#include "Mesh\Vertex.h"
#include "Mesh\Edge.h"
#include "Mesh\Face.h"
#include "Mesh\HalfEdge.h"
#include "Mesh\BaseMesh.h"

#include "Mesh\boundary.h"
#include "Mesh\iterators.h"
#include "Parser\parser.h"

#include "bmp/RgbImage.h"
#include "math.h"

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

namespace MeshLib
{
    class CMyVertex;
    class CMyEdge;
    class CMyFace;
    class CMyHalfEdge;

   
    /*! \brief CMyVertex class
     *
     *   Vertex class for this demo
     *   Trait : Vertex color
     */
    class CMyVertex : public CVertex
    {
    public:
        /*! constructor */
        CMyVertex() : m_rgb(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0) {};
        
        /*! read vertex attributes */
        void _from_string();

        /*! write vertex attributes */
        void _to_string();

        /*! vertex color */
        CPoint& rgb() { return m_rgb; };

        CPoint& geo_rgb() { return m_geo_rgb; };
    protected:
        /*! vertex color */
        CPoint m_rgb;

        /*! vretex xyz -> rgb */
        CPoint m_geo_rgb;
    };

    inline void CMyVertex::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken* token = *iter;
            if (token->m_key == "uv") //CPoint2
            {
                token->m_value >> m_uv;
            }
            if (token->m_key == "rgb") // CPoint
            {
                token->m_value >> m_rgb;
            }
        }
    }

    inline void CMyVertex::_to_string()
    {
        CParser parser(m_string);
        parser._removeToken("uv");

        parser._toString(m_string);
        std::stringstream iss;

        iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

        if (m_string.length() > 0)
        {
            m_string += " ";
        }
        m_string += iss.str();
    }
    
    /*! \brief CMyEdge class
     *
     *   Edge class for this demo
     *   Trait : Edge sharp
     */
    class CMyEdge : public CEdge
    {
    public:
        /*! constructor */
        CMyEdge() :m_sharp(false),m_weight(0) {};

        /*! read edge attributes */
        void _from_string();

        /*! write edge attributes */
        void _to_string();

        /*! sharp edge */
        bool& sharp() { return m_sharp; };

        double& weight() { return m_weight; }

        double& length() { return m_length; }
    protected:
        /*! sharp edge */
        bool m_sharp;
        double m_weight;
        double m_length;
    };

    inline void CMyEdge::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken* token = *iter;
            if (token->m_key == "sharp") // bool
            {
                m_sharp = true;
            }
        }
    }

    inline void CMyEdge::_to_string()
    {
        CParser parser(m_string);
        parser._removeToken("sharp");

        parser._toString(m_string);
        std::stringstream iss;

        if (m_sharp)
            iss << "sharp";

        if (m_string.length() > 0)
        {
            m_string += " ";
        }
        m_string += iss.str();
    }

    /*! \brief CMyFace class
     *
     *   Face class for this demo
     *   Trait : Face normal
     */
    class CMyFace : public CFace
    {
    public:
        /*! face normal */
        CPoint& normal() { return m_normal; };
    protected:
        /*! face normal */
        CPoint m_normal;
    };

    /*! \brief CMyHalfEdge class
     *
     *   HalfEdge class for this demo
     */
    class CMyHalfEdge : public CHalfEdge
    {
    protected:
        double m_angle;
    public:
        double& angle() { return m_angle; }
    };

    /*! \brief CMyMesh class
     *
     *	Mesh class for this demo
     *
     */
    template<typename V, typename E, typename F, typename H>
    class MyMesh : public CBaseMesh<V, E, F, H>
    {
    public:
        typedef V V;
        typedef E E;
        typedef F F;
        typedef H H;

        typedef CBoundary<V, E, F, H>					CBoundary;
        typedef CLoop<V, E, F, H>						CLoop;

        typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
        typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
        typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
        typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;

        typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
        typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
        typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
        typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
        typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;

        typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
        typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
        typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;

        /*!
         *  \brief output the mesh information
         */
        void outputMeshInfo();
        
        /*!
         *  \brief a demo to show how to use iterators of MeshLib
         */
        void testIterator();


        void geometricImage();

    protected:
        F* _locate(CPoint2& p);

        CPoint _interpolation_xyz(CPoint2& p, F* f);
    };

    typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

    template<typename V, typename E, typename F, typename H>
    void MeshLib::MyMesh<V, E, F, H>::outputMeshInfo()
    {
        int nv = this->numVertices();
        int ne = this->numEdges();
        int nf = this->numFaces();

        std::cout << "#V=" << nv << "  ";
        std::cout << "#E=" << ne << "  ";
        std::cout << "#F=" << nf << "  ";

        int euler_char = nv - ne + nf;
        std::cout << "Euler's characteristic=" << euler_char << "  ";

        CBoundary boundary(this);
        std::vector<CLoop*>& loops = boundary.loops();
        int nb = loops.size();

        int genus = (2 - (euler_char + nb)) / 2;
        std::cout << "genus=" << genus << std::endl;
    }

    template<typename V, typename E, typename F, typename H>
    void MyMesh<V, E, F, H>::testIterator()
    {
        for (MeshVertexIterator viter(this); !viter.end(); ++viter)
        {
            V* pV = *viter;
            // you can do something to the vertex here
            // ...

            for (VertexVertexIterator vviter(pV); !vviter.end(); ++vviter)
            {
                     V* pW = *vviter;
                     // you can do something to the neighboring vertices with CCW
                     // ...
             }

            for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
            {
             E* pE = *veiter;
                // you can do something to the neighboring edges with CCW
                // ...
            }

            for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
            {
                F* pF = *vfiter;
                // you can do something to the neighboring faces with CCW
                // ...
            }

            for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
            {
                H* pH = *vhiter;
                // you can do something to the incoming halfedges with CCW
        // ...
}
        }

        for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
        {
            E* pE = *eiter;
            // you can do something to the edge here
            // ...
        }

        for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
        {
            F* pF = *fiter;
            // you can do something to the face here
            // ...
        }

        //there are some other iterators which you can find them in class MyMesh

        std::cout << "Iterators test OK.\n";
    }

    template<typename V, typename E, typename F, typename H>
    void MyMesh<V, E, F, H>::geometricImage()
    {
        const CPoint s(1, 1, 1);
        const CPoint2 s1(1, 1);
        /*
        for (int i = 0; i < 4; i++) {
            cout << this->idVertex(cornerID[0]).uv()[0];
        }
        */
        for (MeshVertexIterator viter(this); !viter.end(); ++viter) {
            V* pV = *viter;
            CPoint& p = pV->point();
            p += s;
            p /= 2;
        }for (MeshVertexIterator viter(this); !viter.end(); ++viter) {
            V* pV = *viter;
            CPoint2& p = pV->uv();
            p += s1;
            p /= 2;
        }
        
        const int w = 640;
        const int h = 640;

        RgbImage image(h, w);

        CPoint2 a(0, 0);
        const double width = 1.0;
        const double height = 1.0;

        double step = width / (w - 1);
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                CPoint2 p(0, 0);
                p[0] = a[0] + j * step;
                p[1] = a[1] + i * step;

                F* pF = _locate(p);
                if (pF != NULL) {
                    double t, k;
                    CPoint xyz = _interpolation_xyz(p, pF);
                    image.SetRgbPixelf(i, j, xyz[0], xyz[1], xyz[2]);
                }
                else {
                    printf("error: face not found");
                }
            }
            
        }
        image.WriteBmpFile("geometric_image.bmp");
    }
    
    template<typename V, typename E, typename F, typename H>
    F* MyMesh<V, E, F, H>::_locate(CPoint2& p)
    {
        F* pF = NULL;

        for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter) {
            pF = *fiter;
            V* trivertex[3];
            double area[4];
            int i = 0;
            double edge_len[6];//AB,BC,AC,pA,pB,pC
            for (FaceVertexIterator fviter(pF); !fviter.end(); ++fviter) {
                trivertex[i] = *fviter;
                i++;
            }

            edge_len[0] = (trivertex[0]->uv() - trivertex[1]->uv()).norm();//AB
            edge_len[1] = (trivertex[1]->uv() - trivertex[2]->uv()).norm();//BC
            edge_len[2] = (trivertex[0]->uv() - trivertex[2]->uv()).norm();//AC
            edge_len[3] = (trivertex[0]->uv() - p).norm();//pA
            edge_len[4] = (trivertex[0]->uv() - p).norm();//pB
            edge_len[5] = (trivertex[0]->uv() - p).norm();//pC

            area[0] = Sarea(edge_len[0], edge_len[1], edge_len[2]);
            area[1] = Sarea(edge_len[3], edge_len[4], edge_len[0]);
            area[2] = Sarea(edge_len[3], edge_len[5], edge_len[2]);
            area[3] = Sarea(edge_len[4], edge_len[5], edge_len[1]);

            if (area[0] - area[1] - area[2] - area[3] < 0.001) {
                return pF;
            }
            
        }               
        return NULL;
    }

    template<typename V, typename E, typename F, typename H>
    CPoint MyMesh<V, E, F, H>::_interpolation_xyz(CPoint2& p, F* f)
    {
        H* pH = this->faceHalfedge(f);
        V* pA = this->halfedgeSource(pH);
        V* pB = this->halfedgeTarget(pH);
        pH = this->halfedgeNext(pH);
        V* pC = this->halfedgeTarget(pH);
        
        double t = 0, k = 0;
        t = ((p[0] - pC->uv()[0]) * (pB->uv()[1] - pC->uv()[1]) - (pB->uv()[0] - pC->uv()[0]) * (p[1] - pC->uv()[1])) / ((pA->uv()[0] - pC->uv()[0]) * (pB->uv()[1] - pC->uv()[1]) - (pB->uv()[0] - pC->uv()[0]) * (pA->uv()[1] - pC->uv()[1]));
        k = ((pA->uv()[0] - pC->uv()[0]) * (p[1] - pC->uv()[1]) - (p[0] - pC->uv()[0]) * (pA->uv()[1] - pC->uv()[1])) / ((pA->uv()[0] - pC->uv()[0]) * (pB->uv()[1] - pC->uv()[1]) - (pB->uv()[0] - pC->uv()[0]) * (pA->uv()[1] - pC->uv()[1]));
       
        CPoint xyz = pA->point() * t + pB->point() * k + pC->point() * (1 - t - k);
        return xyz;
    }
    double Sarea(double a, double b, double c) {
        double s = (a + b + c) / 2;
        return sqrt(s * (s - a) * (s - b) * (s - c));
    }

}

template <typename M>
class CHarmonicMapper {
public:
    void _set_boundary();
    void _iterative_map(double epsilon);
    CHarmonicMapper(M* m) :m_pMesh(m) {
        _calculate_edge_weight();
    }
    void _calculate_edge_weight();
    double _inverse_cosine_law(double a, double b, double c) {
        return acos((b * b + c * c - a * a) / (2 * b * c));
    }
protected:
    M* m_pMesh;
    MeshLib::CMyVertex* corner[4];
};

template <typename M>
void CHarmonicMapper<M>::_iterative_map(double epsilon) {
    _set_boundary();

    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        MeshLib::CMyVertex* pV = *viter;
        if (pV->boundary()) { 
            //std::cout << pV->uv()[0] <<" "<< pV->uv()[1] << std::endl;
            continue; 
        }
        pV->uv() = CPoint2(0, 0);
    }
    while (true) {
        double error = -1e+10;

        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {

            MeshLib::CMyVertex* pV = *viter;
            if (pV->boundary()) continue;

            double sw = 0;
            CPoint2 suv(0, 0);
            for (M::VertexVertexIterator vviter(pV); !vviter.end(); ++vviter) {
                MeshLib::CMyVertex* pW = *vviter;
                MeshLib::CMyEdge* pE = m_pMesh->vertexEdge(pV, pW);
                double w = pE->weight();
                //std::cout << w << std::endl;
                sw += w;
                suv = suv + pW->uv() * w;
            }
            suv /= sw;
            double verror = (pV->uv() - suv).norm();
            //std::cout << suv[0]<<" "<<suv[1] << std::endl;
            error = (verror > error) ? verror : error;
            pV->uv() = suv;
        }
        printf("Current max error is %f\n", error);
        if (error < epsilon) break;
    }
}
template <typename M>
void CHarmonicMapper<M>::_calculate_edge_weight() {
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
        CMyEdge* pE = *eiter;
        CMyVertex* v1 = m_pMesh->edgeVertex1(pE);
        CMyVertex* v2 = m_pMesh->edgeVertex2(pE);
        pE->length() = (v1->point() - v2->point()).norm();
    }

    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter) {
        CMyFace* pF = *fiter;
        CMyHalfEdge* pH[3];
        pH[0] = m_pMesh->faceHalfedge(pF);

        for (int i = 0; i < 3; i++) {
            pH[(i + 1) % 3] = m_pMesh->faceNextCcwHalfEdge(pH[i]);
        }
        double len[3];
        for (int i = 0; i < 3; i++) {
            len[i] = m_pMesh->halfedgeEdge(pH[i])->length();
        }
        for (int i = 0; i < 3; i++) {
            double a = len[(i + 1) % 3], b = len[(i + 2) % 3], c = len[i];
            pH[(i + 1) % 3]->angle() = _inverse_cosine_law(a, b, c);
        }
    }

    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
        CMyEdge* pE = *eiter;
        if (!pE->boundary()) {
            double theta[2];
            theta[0] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            theta[1] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 1))->angle();
            pE->weight() = cos(theta[0]) / sin(theta[0]) + cos(theta[1]) / sin(theta[1]);
        }
        else {
            double theta = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            pE->weight() = cos(theta) / sin(theta);
        }
    }
}
template<typename M>
void CHarmonicMapper<M>::_set_boundary() {
    double len[4] = { 0, 0, 0, 0 };
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        CMyVertex* pV = *viter;
        if (!m_pMesh->isBoundary(pV)) continue;
        if (pV->point()[0] > 0 && pV->point()[1] > 0 && pV->point().norm() > len[0]) {
            len[0] = pV->point().norm();
            corner[0] = pV;
        }
        if (pV->point()[0] > 0 && pV->point()[1] < 0 && pV->point().norm() > len[1]) {
            len[1] = pV->point().norm();
            corner[1] = pV;
        }
        if (pV->point()[0] < 0 && pV->point()[1] < 0 && pV->point().norm()> len[2]) {
            len[2] = pV->point().norm();
            corner[2] = pV;
        }
        if (pV->point()[0] < 0 && pV->point()[1] > 0 && pV->point().norm() > len[3]) {
            len[3] = pV->point().norm();
            corner[3] = pV;
        }
    }
    double len1[4] = { 0, 0, 0, 0 };
    double len2[4] = { 0, 0, 0, 0 };
    for (int i = 0; i < 4; i++) {
        CMyHalfEdge* he = m_pMesh->vertexMostCcwInHalfEdge(corner[i]);
        do {
            CMyVertex* v = (CMyVertex*)he->target();
            he = m_pMesh->vertexMostClwOutHalfEdge(v);
            len1[i] += m_pMesh->edgeLength((CMyEdge*)he->edge());
        } while (he->source() != corner[(i + 1) % 4]);
    }
    for (int i = 0; i < 4; i++) {
        CMyHalfEdge* he = m_pMesh->vertexMostCcwInHalfEdge(corner[i]);
        do {
            CMyVertex* v = (CMyVertex*)he->target();
            he = m_pMesh->vertexMostClwOutHalfEdge(v);
            len2[i] += m_pMesh->edgeLength((CMyEdge*)he->edge());
            switch (i) {
            case 0:
                v->uv()[0] = 1;
                v->uv()[1] = 1 - len2[i] / len1[i] * 2;
                break;
            case 1:
                v->uv()[0] = 1 - len2[i] / len1[i] * 2;
                v->uv()[1] = -1;
                break;
            case 2:
                v->uv()[0] = -1;
                v->uv()[1] = -1 + len2[i] / len1[i] * 2;
                break;
            case 3:
                v->uv()[0] = -1 + len2[i] / len1[i] * 2;
                v->uv()[1] = 1;
            }
        } while (he->source() != corner[(i + 1) % 4]);
    }
    corner[0]->uv() = CPoint2(1, 1);
    corner[1]->uv() = CPoint2(1, -1);
    corner[2]->uv() = CPoint2(-1, -1);
    corner[3]->uv() = CPoint2(-1, 1);
}





#endif // !_MY_MESH_
