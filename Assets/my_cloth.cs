using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class my_cloth : MonoBehaviour
{
    //
    public int resolution = 11;

    // Mass of each point in the cloth
    public float pointMass = 0.1f;

    // Spring constant for each connection between points in the cloth
    public float sturctSpringConstant = -50.0f;
    public float shearSpringConstant = 0.0f;
    public float bendSpringConstant = 0.0f;

    // Damping factor for each point in the cloth
    public float dampingConstant = -0.1f;

    // Strength and direction of gravity acting on the cloth
    public Vector3 gravity = new Vector3(0.0f, -9.81f, 0.0f);
    public float drag = 0.0f;

    // Strength and direction of wind affecting the cloth
    public Vector3 wind = Vector3.zero;

    // Whether to enable collision detection with other objects in the scene
    public bool collisionDetection = true;

    // Integration method used to calculate the position and velocity of each point in the cloth over time
    // public IntegrationMethod integrationMethod = IntegrationMethod.Verlet;

    // Mesh data for the cloth object
    private MeshFilter meshFilter;
    private MeshCollider meshCollider;

    private Mesh clothMesh;
    private Mesh originalCopy;
    private Mesh mesh;

    private Vector3[] vertices;
    private Vector3[] verticesCopy;
    private Vector3[] vertexVelocities;
    private float[] vertexMass;
    private Vector3[] nextPos;
    private Vector3[] nextVel;

    private int[] triangles;
    List<Edge> edges = new List<Edge>();
    private Vector3[] forces;

    public struct Edge
    {
        public int vertexIndexA;
        public int vertexIndexB;
        public int flag; // 0: struct, 1: shear, 2: bend

        public Edge(int a, int b, int f)
        {
            vertexIndexA = Mathf.Min(a, b);
            vertexIndexB = Mathf.Max(a, b);
            flag = f;
        }

        public override bool Equals(object obj)
        {
            Edge other = (Edge)obj;
            return (vertexIndexA == other.vertexIndexA && vertexIndexB == other.vertexIndexB);
        }

        public override int GetHashCode()
        {
            return vertexIndexA ^ vertexIndexB;
        }
    }

    // Constraints on the motion of the cloth
    // public Constraint[] constraints;

    // Other variables and methods for simulating cloth physics...
    // Start is called before the first frame update
    void Start()
    {
        mesh = GetComponent<MeshFilter>().mesh;

        // Resize the mesh.
        int n = resolution;
        Vector3[] X = new Vector3[n * n]; // new vertices
        Vector2[] UV = new Vector2[n * n]; // UV
        int[] triangles = new int[(n - 1) * (n - 1) * 6]; // new triangle
        // vertices position setup
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
            {
                X[j * n + i] = new Vector3(5 - 10.0f * i / (n - 1), 0, 5 - 10.0f * j / (n - 1));
                UV[j * n + i] = new Vector3(i / (n - 1.0f), j / (n - 1.0f));
            }
        // triangle setup
        int t = 0;
        for (int j = 0; j < n - 1; j++)
            for (int i = 0; i < n - 1; i++)
            {
                triangles[t * 6 + 5] = j * n + i;
                triangles[t * 6 + 4] = j * n + i + 1;
                triangles[t * 6 + 3] = (j + 1) * n + i + 1;
                triangles[t * 6 + 2] = j * n + i;
                triangles[t * 6 + 1] = (j + 1) * n + i + 1;
                triangles[t * 6 + 0] = (j + 1) * n + i;
                t++;
            }
        //  re-define
        mesh.vertices = X;
        mesh.triangles = triangles;
        mesh.uv = UV;
        mesh.RecalculateNormals();

        meshFilter = GetComponent<MeshFilter>();
        meshCollider = GetComponent<MeshCollider>();

        clothMesh = meshFilter.mesh;

        Mesh clothCopy = new Mesh();
        clothCopy.name = "Cloth";
        clothCopy.vertices = clothMesh.vertices;
        clothCopy.triangles = clothMesh.triangles;
        clothCopy.uv = clothMesh.uv;
        clothCopy.RecalculateNormals();
        GetComponent<MeshFilter>().mesh = clothCopy;
        GetComponent<MeshCollider>().sharedMesh = clothCopy;

        // Check if the mesh filter component exists
        if (meshFilter != null)
        {
            // Mesh
            clothMesh = GetComponent<MeshFilter>().mesh;
            originalCopy = meshFilter.mesh;
            // Mesh data
            vertices = clothMesh.vertices;
            triangles = clothMesh.triangles;

            vertexVelocities = new Vector3[vertices.Length];
            vertexMass = new float[vertices.Length];
            forces = new Vector3[vertices.Length];
            nextPos = new Vector3[vertices.Length];
            nextVel = new Vector3[vertices.Length];

            // original copy
            verticesCopy = originalCopy.vertices;

            // init
            for (int i = 0; i < vertices.Length; i++)
            {
                vertexVelocities[i] = Vector3.zero;
                vertexMass[i] = pointMass;
                forces[i] = Vector3.zero;
                nextPos[i] = vertices[i];
                nextVel[i] = Vector3.zero;
            }

            // unique edge init
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    // [structual] right edge
                    if (j < resolution - 1)
                    {
                        Edge edge1 = new Edge(i * resolution + j, i * resolution + j + 1, 0);
                        AddEdgeIfUnique(edges, edge1);
                    }
                    // [structual] down edge
                    if (i < resolution - 1)
                    {
                        Edge edge2 = new Edge(
                            i * resolution + j,
                            i * resolution + j + resolution,
                            0
                        );
                        AddEdgeIfUnique(edges, edge2);
                    }
                    // [shearing] right bottom edge
                    if (j < resolution - 1 && i < resolution - 1)
                    {
                        Edge edge3 = new Edge(
                            i * resolution + j,
                            i * resolution + j + resolution + 1,
                            1
                        );
                        // Debug.Log(i * resolution+j );
                        // Debug.Log( i * resolution + j + resolution + 1);
                        AddEdgeIfUnique(edges, edge3);
                    }
                    // [shearing] left bottom edge
                    if (j % resolution != 0 && i < resolution - 1)
                    {
                        Edge edge4 = new Edge(
                            i * resolution + j,
                            i * resolution + j + resolution - 1,
                            1
                        );
                        AddEdgeIfUnique(edges, edge4);
                    }
                    // [bending] right edge
                    if (j < resolution - 2)
                    {
                        Edge edge5 = new Edge(i + j, i + j + 2, 2);
                        AddEdgeIfUnique(edges, edge5);
                    }
                    // [bending] down edge
                    if (i < resolution - 2)
                    {
                        Edge edge6 = new Edge(i + j, i + j + resolution + 2, 2);
                        AddEdgeIfUnique(edges, edge6);
                    }
                }
            }
        }
    }

    // Update is called once per frame
    void FixedUpdate()
    // void Update()
    {
        // reset forces
        for (int i = 0; i < vertices.Length; i++)
        {
            forces[i] = Vector3.zero;
        }

        // TODO:
        // External Force (Gravity, Drag)
        // Internal Force (Spring, Damping)
        // Vertices Update (Position, Velocity)

        ExternalForce();
        InternalForce();

        // Collision detection
        if (collisionDetection)
        {
            CollisionDetection();
        }

        VerticesUpdate();
        // Update the mesh vertices
    }

    void AddEdgeIfUnique(List<Edge> edges, Edge edge)
    {
        bool isUnique = true;

        foreach (Edge e in edges)
        {
            if (e.Equals(edge))
            {
                isUnique = false;
                break;
            }
        }

        if (isUnique)
        {
            edges.Add(edge);
        }
    }

    private void ExternalForce()
    {
        for (int i = 0; i < vertices.Length; i++)
        {
            // F = mass * accerlaration (gravity)
            // F = drag * velocity
            forces[i] += vertexMass[i] * gravity;
            forces[i] += vertexVelocities[i] * drag;
        }
    }

    private void InternalForce()
    {
        // thru unique edges instead of triangle

        foreach (Edge e in edges)
        {
            int v1 = e.vertexIndexA;
            int v2 = e.vertexIndexB;
            int f = e.flag;

            forces[v1] += CalculateSpringAndDamperForce(v1, v2, f);
            forces[v2] += CalculateSpringAndDamperForce(v2, v1, f);
        }
    }

    private Vector3 CalculateSpringAndDamperForce(int vi1, int vi2, int flag)
    {
        float springConstant = 0;
        switch (flag)
        {
            case 0:
                springConstant = sturctSpringConstant;
                break;
            case 1:
                springConstant = shearSpringConstant;
                break;
            case 2:
                springConstant = bendSpringConstant;
                break;
            default:
                break;
        }

        // Get the current positions of the two vertices
        Vector3 vertex1Position = vertices[vi1];
        Vector3 vertex2Position = vertices[vi2];

        // Calculate the vector between the two vertices
        Vector3 dir = vertex2Position - vertex1Position;

        // Calculate the distance between the two vertices
        float length = dir.magnitude;

        // Calculate the direction of the spring force
        Vector3 normalizedDir = dir.normalized;

        // Get the restLength from copy
        Vector3 originalDir = verticesCopy[vi2] - verticesCopy[vi1];
        float restLength = originalDir.magnitude;

        // Calculate the spring force using Hooke's Law (F = -kx)
        float springForceMagnitude = -springConstant * (length - restLength);
        Vector3 springForce = springForceMagnitude * normalizedDir;

        // Calculate the damper force using a damping coefficient (F = -bv)
        Vector3 relativeVelocity = vertexVelocities[vi2] - vertexVelocities[vi1];
        float damperForceMagnitude =
            -dampingConstant * Vector3.Dot(relativeVelocity, normalizedDir);
        Vector3 damperForce = damperForceMagnitude * normalizedDir;

        return springForce + damperForce;
    }

    private void VerticesUpdate()
    {
        //
        for (int i = 0; i < vertices.Length; i++)
        {
            if (i != 0 && i != resolution - 1)
            {
                Vector3 acc = forces[i] / vertexMass[i];
                vertexVelocities[i] += acc * Time.deltaTime;
                vertices[i] += vertexVelocities[i] * Time.deltaTime;
            }
        }
        clothMesh.vertices = vertices;
        clothMesh.RecalculateNormals();

        meshCollider.sharedMesh = clothMesh;
    }

    private void OnCollisionEnter(Collision collision)
    {
        ContactPoint[] contacts = collision.contacts;
        // Debug.Log("Number of contacts: " + contacts.Length);

        foreach (ContactPoint contact in contacts)
        {
            // Debug.Log("Collided with " + contact.otherCollider.gameObject.name);
            // Mesh otherMesh = contact.otherCollider.GetComponent<MeshFilter>().mesh;
            // Vector3 p = contact.otherCollider.transform.InverseTransformPoint(contact.point);
            // int triangle = ClosestPointOnMesh(otherMesh, p);
            // Debug.Log("Collided with " + triangle);
            // for (int i = 0; i < otherMesh.triangles.Length; i += 3)
            // {
            //     int index1 = otherMesh.triangles[i];
            //     int index2 = otherMesh.triangles[i + 1];
            //     int index3 = otherMesh.triangles[i + 2];

            //     Vector3 vertex1 = otherMesh.vertices[index1];
            //     Vector3 vertex2 = otherMesh.vertices[index2];
            //     Vector3 vertex3 = otherMesh.vertices[index3];

            //     Vector3 p = contact.otherCollider.transform.InverseTransformPoint(contact.point);

            //     Vector3 wv1 = contact.otherCollider.transform.TransformPoint(vertex1);
            //     Vector3 wv2 = contact.otherCollider.transform.TransformPoint(vertex2);
            //     Vector3 wv3 = contact.otherCollider.transform.TransformPoint(vertex3);

            //     if (PointInTriangle(p, vertex1, vertex2, vertex3))
            //     {
            //         // Vector3 side1 = vertex1 - vertex2;
            //         // Vector3 side2 = vertex3 - vertex2;
            //         // Vector3 normal = Vector3.Cross(side1, side2).normalized;
            //         // Debug.Log("normal is: " + normal);
            //         Debug.Log("Triangle: " + i / 3);
            //     }
            // }
        }
    }

    int ClosestPointOnMesh(Mesh mesh, Vector3 point)
    {
        float closestDistance = Mathf.Infinity;
        Vector3 closestPoint = Vector3.zero;
        int tri = -1;

        // Loop through each triangle of the mesh
        for (int i = 0; i < mesh.triangles.Length; i += 3)
        {
            // Get the vertices of the triangle
            Vector3 v1 = mesh.vertices[mesh.triangles[i]];
            Vector3 v2 = mesh.vertices[mesh.triangles[i + 1]];
            Vector3 v3 = mesh.vertices[mesh.triangles[i + 2]];

            // Find the closest point on the triangle to the input point
            float distance = DistanceToTriangle(point, v1, v2, v3);

            // If the distance is smaller than the current closest distance, update the closest point
            if (distance < closestDistance)
            {
                closestDistance = distance;
                // closestPoint = triangleClosestPoint;
                tri = i / 3;
            }
        }

        return tri;
    }

    // Returns the closest point on a triangle to a given point
    float DistanceToTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
    {
        Vector3 normal = Vector3.Cross(b - a, c - a).normalized;
        float distance = Vector3.Dot(normal, a - p);
        Vector3 projectedPoint = p + distance * normal;
        if (PointInTriangle(projectedPoint, a, b, c))
        {
            return Vector3.Distance(projectedPoint, p);
        }
        else
        {
            float minDistance = Mathf.Infinity;
            minDistance = Mathf.Min(
                minDistance,
                Vector3.Distance(p, ClosestPointOnSegment(a, b, projectedPoint))
            );
            minDistance = Mathf.Min(
                minDistance,
                Vector3.Distance(p, ClosestPointOnSegment(b, c, projectedPoint))
            );
            minDistance = Mathf.Min(
                minDistance,
                Vector3.Distance(p, ClosestPointOnSegment(c, a, projectedPoint))
            );
            return minDistance;
        }
    }

    public static Vector3 ClosestPointOnSegment(Vector3 a, Vector3 b, Vector3 p)
    {
        Vector3 ab = b - a;
        float t = Vector3.Dot(p - a, ab) / Vector3.Dot(ab, ab);
        t = Mathf.Clamp01(t);
        return a + t * ab;
    }

    private bool PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
    {
        // Check if point is inside the triangle defined by a, b, and c
        Vector3 v0 = c - a;
        Vector3 v1 = b - a;
        Vector3 v2 = p - a;

        float dot00 = Vector3.Dot(v0, v0);
        float dot01 = Vector3.Dot(v0, v1);
        float dot02 = Vector3.Dot(v0, v2);
        float dot11 = Vector3.Dot(v1, v1);
        float dot12 = Vector3.Dot(v1, v2);

        float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= 0) && (v >= 0) && (u + v < 1);
    }

    private void CollisionDetection()
    {
        MeshFilter[] meshFilters = FindObjectsOfType<MeshFilter>();

        // for (int i = 0; i < clothMesh.triangles.Length; i += 3)
        // // for (int i = 0; i < 1; i += 3)
        // {
        //     int index1 = clothMesh.triangles[i];
        //     int index2 = clothMesh.triangles[i + 1];
        //     int index3 = clothMesh.triangles[i + 2];

        //     Vector3 vertex1 = clothMesh.vertices[index1];
        //     Vector3 vertex2 = clothMesh.vertices[index2];
        //     Vector3 vertex3 = clothMesh.vertices[index3];
        //     //  in world coord
        //     Vector3 wv1 = meshFilter.transform.TransformPoint(vertex1);
        //     Vector3 wv2 = meshFilter.transform.TransformPoint(vertex2);
        //     Vector3 wv3 = meshFilter.transform.TransformPoint(vertex3);

        //     foreach (MeshFilter filter in meshFilters)
        //     {
        //         if (filter.transform != transform)
        //         {
        //             Debug.Log(filter.gameObject.name);
        //             Mesh other = filter.mesh;
        //             for (int j = 0; j < other.triangles.Length; j += 3)
        //             {
        //                 int OtherIndex1 = other.triangles[j];
        //                 int OtherIndex2 = other.triangles[j + 1];
        //                 int OtherIndex3 = other.triangles[j + 2];

        //                 Vector3 Othervertex1 = other.vertices[OtherIndex1];
        //                 Vector3 Othervertex2 = other.vertices[OtherIndex2];
        //                 Vector3 Othervertex3 = other.vertices[OtherIndex3];
        //                 //  in world coord
        //                 Vector3 Otherwv1 = filter.transform.TransformPoint(Othervertex1);
        //                 Vector3 Otherwv2 = filter.transform.TransformPoint(Othervertex2);
        //                 Vector3 Otherwv3 = filter.transform.TransformPoint(Othervertex3);

        //                 if (IsTrianglesIntersect(wv1, wv2, wv3, Otherwv1, Otherwv2, Otherwv3))
        //                 {
        //                     Debug.Log("collided");
        //                     Plane plane = new Plane(Otherwv1, Otherwv2, Otherwv3);
        //                     Vector3 hit_normal = plane.normal;
        //                     Vector3 vN, vT, v, temp_force;

        //                     // normal force from the collision object
        //                     vN = Vector3.Dot(nextVel[index1], hit_normal) * hit_normal;
        //                     vT = nextVel[index1] - vN;
        //                     v = -0.001f * vN + vT;
        //                     // recalculate force
        //                     temp_force =
        //                         (vertexMass[index1] * (v - vertexVelocities[index1])) / Time.deltaTime;
        //                     forces[index1] = temp_force;

        //                     // normal force from the collision object
        //                     vN = Vector3.Dot(nextVel[index2], hit_normal) * hit_normal;
        //                     vT = nextVel[index2] - vN;
        //                     v = -0.001f * vN + vT;
        //                     // recalculate force
        //                     temp_force =
        //                         (vertexMass[index2] * (v - vertexVelocities[index2])) / Time.deltaTime;
        //                     forces[index2] = temp_force;

        //                     // normal force from the collision object
        //                     vN = Vector3.Dot(nextVel[index3], hit_normal) * hit_normal;
        //                     vT = nextVel[index3] - vN;
        //                     v = -0.001f * vN + vT;
        //                     // recalculate force
        //                     temp_force =
        //                         (vertexMass[index3] * (v - vertexVelocities[index3])) / Time.deltaTime;
        //                     forces[index3] = temp_force;

        //                 }
        //             }
        //         }

        //         // Debug.Log(filter.mesh.name);
        //     }
        // }

        float thickness = 1e-5f;
        // particle collision (it will interpenatrate because of edges)
        float rayDistance = 0;
        Vector3 rayDirection;
        // RaycastHit hit;
        for (int i = 0; i < vertices.Length; i++)
        {
            Vector3 acc = forces[i] / vertexMass[i];
            nextVel[i] = vertexVelocities[i] + acc * Time.deltaTime;
            nextPos[i] = vertices[i] + nextVel[i] * Time.deltaTime;

            rayDistance = Vector3.Distance(vertices[i], nextPos[i]);
            rayDirection = (nextPos[i] - vertices[i]).normalized;

            // particle collision
            RaycastHit[] hits = Physics.RaycastAll(
                meshFilter.transform.TransformPoint(vertices[i] - thickness * rayDirection),
                rayDirection,
                rayDistance,
                -1
            );

            foreach (RaycastHit hit in hits)
            {
                if (hit.transform != transform)
                {
                    // normal force from the collision object
                    Vector3 vN = Vector3.Dot(nextVel[i], hit.normal) * hit.normal;
                    Vector3 vT = nextVel[i] - vN;
                    Vector3 v = -0.001f * vN + vT;
                    // recalculate force
                    // Vector3 temp_force =
                    //     (vertexMass[i] * (v - vertexVelocities[i])) / Time.deltaTime;
                    // forces[i] = temp_force;
                    break;
                }
            }
        }

        foreach (Edge e in edges)
        {
            Edgecollider(e.vertexIndexA, e.vertexIndexB);
            // if(e.flag == 0)
            // {
            //     // Edgecollider(e.vertexIndexA, e.vertexIndexB);
            // }

        }
    }

    private bool IsTrianglesIntersect(
        Vector3 p1,
        Vector3 p2,
        Vector3 p3,
        Vector3 q1,
        Vector3 q2,
        Vector3 q3
    )
    {
        // Compute the planes containing the triangles
        Plane plane1 = new Plane(p1, p2, p3);
        Plane plane2 = new Plane(q1, q2, q3);
        Vector3 n1 = plane1.normal;
        Vector3 n2 = plane2.normal;
        float dot = Vector3.Dot(n1, n2);

        // Check if the planes are parallel
        if (Mathf.Abs(dot - 1.0f) < 0.0001f)
        {
            Debug.Log("parallel");
            return false;
        }

        // Compute the intersection line between the planes
        Vector3 lineDirection = Vector3.Cross(plane1.normal, plane2.normal);
        Vector3 linePoint = Vector3.zero;
        if (!LinePlaneIntersection(plane1, plane2, out linePoint, out lineDirection))
        {
            Debug.Log("not intersect");
            return false;
        }

        // Check if the line intersects the triangles' edges
        if (
            LineTriangleIntersection(linePoint, lineDirection, p1, p2, p3)
            || LineTriangleIntersection(linePoint, lineDirection, q1, q2, q3)
        )
        {
            Debug.Log("impossible1");
            return true;
        }

        // Check if the line intersects the triangles' planes
        if (PointInTriangle(linePoint, p1, p2, p3) && PointInTriangle(linePoint, q1, q2, q3))
        {
            Debug.Log("impossible2");
            return true;
        }

        Debug.Log("f");
        return false;
    }

    private bool LinePlaneIntersection(
        Plane plane1,
        Plane plane2,
        out Vector3 linePoint,
        out Vector3 lineDirection
    )
    {
        lineDirection = Vector3.Cross(plane1.normal, plane2.normal);
        float denom = Vector3.Dot(lineDirection, lineDirection);

        if (Mathf.Approximately(denom, 0))
        {
            linePoint = Vector3.zero;
            return false;
        }

        Vector3 diff = plane2.normal * plane1.distance - plane1.normal * plane2.distance;
        linePoint = Vector3.Cross(diff, lineDirection) / denom;
        return true;
    }

    private bool LineTriangleIntersection(
        Vector3 linePoint,
        Vector3 lineDirection,
        Vector3 p1,
        Vector3 p2,
        Vector3 p3
    )
    {
        if (PointOnTriangleEdges(linePoint, p1, p2, p3))
        {
            return true;
        }

        // Check if the line intersects the triangle plane
        Plane trianglePlane = new Plane(p1, p2, p3);
        float t =
            (trianglePlane.distance - Vector3.Dot(trianglePlane.normal, linePoint))
            / Vector3.Dot(trianglePlane.normal, lineDirection);

        if (t < 0 || t > 1)
        {
            return false;
        }

        // Check if the intersection point is inside the triangle
        Vector3 intersectionPoint = linePoint + t * lineDirection;
        return PointInTriangle(intersectionPoint, p1, p2, p3);
    }

    private bool PointOnTriangleEdges(Vector3 point, Vector3 p1, Vector3 p2, Vector3 p3)
    {
        return PointOnLine(point, p1, p2)
            || PointOnLine(point, p1, p3)
            || PointOnLine(point, p2, p3);
    }

    private bool PointOnLine(Vector3 point, Vector3 lineStart, Vector3 lineEnd)
    {
        float dist = Vector3.Distance(
            point,
            Vector3.Lerp(
                lineStart,
                lineEnd,
                Vector3.Dot(point - lineStart, lineEnd - lineStart)
                    / (lineEnd - lineStart).sqrMagnitude
            )
        );
        return dist < 0.0001;
    }

    private void Edgecollider(int start, int end)
    {
        RaycastHit hit1,
            hit2;

        Vector3 acc1 = forces[start] / vertexMass[start];
        Vector3 acc2 = forces[end] / vertexMass[end];

        Vector3 startVel = vertexVelocities[start] + acc1 * Time.deltaTime;
        Vector3 endVel = vertexVelocities[end] + acc2 * Time.deltaTime;

        Vector3 startPos = vertices[start] + startVel * Time.deltaTime;
        Vector3 endPos = vertices[end] + endVel * Time.deltaTime;

        Vector3 a = meshFilter.transform.TransformPoint(startPos);
        Vector3 b = meshFilter.transform.TransformPoint(endPos);

        float rayDistance = Vector3.Distance(a, b);
        Vector3 rayDirection = (b - a).normalized;
        Vector3 temp_normal = Vector3.zero;

        // particle collision
        RaycastHit[] hits1 = Physics.RaycastAll(a, rayDirection, rayDistance, -1);
        RaycastHit[] hits2 = Physics.RaycastAll(b, -rayDirection, rayDistance, -1);
        bool triggered = false;

        foreach (RaycastHit hit in hits1)
        {
            if (hit.transform != transform)
            {
                temp_normal += hit.normal;
                triggered = true;
                break;
            }
        }
        foreach (RaycastHit hit in hits2)
        {
            if (hit.transform != transform)
            {
                temp_normal += hit.normal;
                triggered = true;
                break;
            }
        }

        temp_normal = temp_normal.normalized;
        Vector3 vN,
            vT,
            v,
            temp_force;

        if (triggered)
        {
            // normal force from the collision object
            vN = Vector3.Dot(nextVel[start], temp_normal) * temp_normal;
            vT = nextVel[start] - vN;
            v = -0.001f * vN + vT;
            // recalculate force
            temp_force = (vertexMass[start] * (v - vertexVelocities[start])) / Time.deltaTime;
            forces[start] = temp_force;

            // normal force from the collision object
            vN = Vector3.Dot(nextVel[end], temp_normal) * temp_normal;
            vT = nextVel[end] - vN;
            v = -0.001f * vN + vT;
            // recalculate force
            temp_force = (vertexMass[end] * (v - vertexVelocities[end])) / Time.deltaTime;
            forces[end] = temp_force;
        }
    }
}
