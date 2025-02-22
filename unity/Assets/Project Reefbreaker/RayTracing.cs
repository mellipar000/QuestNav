using System;
using System.Numerics;
using NetworkTables;
using Unity.Mathematics;
using UnityEngine;
using Matrix4x4 = UnityEngine.Matrix4x4;
using Vector3 = UnityEngine.Vector3;
using Vector4 = UnityEngine.Vector4;

public class RayTracing : MonoBehaviour
{
    private Nt4Client _ntClient;
    private double robotX = 0.0, robotY = 0.0, robotZ = 0.0, robotRotation = 0.0;
    private int[] coralCounters = new int[48], algaeCounters = new int[6], emptyCoralCounters = new int[48], emptyAlgaeCounters = new int[6];
    private bool[] coralPresence = new bool[12], algaePresence = new bool[12], hasSeenCoral = new bool[12], hasSeenAlgae = new bool[12], coralVisibilities = new bool[0], algaeVisibilities = new bool[0]; // presence arrays iterate through all heights of each XY before moving on to the next XY
    private double[] coralDistanceEstimates = new double[0], algaeDistanceEstimates = new double[0];
    private double[][] cameraOffsets = new double[0][], cameraPositions = new double[0][], coralValidIntersections = new double[0][], algaeValidIntersections = new double[0][], coralCorrectIntersections = new double[0][], algaeCorrectIntersections = new double[0][];
    private int[][] cameraParameters = new int[0][];
    private double[][][] coralCameraAngles = new double[0][][], algaeCameraAngles = new double[0][][], coralRayConstants = new double[0][][], algaeRayConstants = new double[0][][], coralRawIntersections = new double[0][][], algaeRawIntersections = new double[0][][];
    private String alliance = "Blue";
    private int _cameraOffsetsSubscriptionId, _coralCameraAnglesSubscriptionId, _algaeCameraAnglesSubscriptionId, _robotPositionSubscriptionId, _allianceSubscriptionId, _cameraParametersSubscriptionId, framePresenceThreshold = 40;

    // Field Constants
    private const double blueReefX = 4.48945, reedReefX = 13.065, reefY = 4.0259, radius = 0.658;
    private const double L1 = 0.4572, L2 = 0.70358, L3 = 1.106678, L4 = 1.652778;
    private readonly double[] coralHeights = new double[] {L1, L2, L3, L4};
    private readonly double[,] blueCoral = new double[,] {
        {3.822, 4.2}, {3.822, 3.859}, // A & B
        {4.0, 3.547}, {4.319, 3.372}, // C & D
        {4.66, 3.372}, {4.963, 3.557}, // E & F
        {5.138, 3.859}, {5.138, 4.2}, // G & H
        {4.963, 4.522}, {4.66, 4.678}, // I & J
        {4.319, 4.678}, {4.0, 4.522} // K & L
    };
    private readonly double[,] redCoral = new double[,] {
        {13.718, 3.859}, {13.718, 4.2}, // A & B
        {13.54, 4.522}, {13.221, 4.678}, // C & D
        {12.88, 4.678}, {12.577, 4.522}, // E & F
        {12.402, 4.2}, {12.402, 3.859}, // G & H
        {12.577, 3.557}, {12.88, 3.372}, // I & J
        {13.221, 3.372}, {13.54, 3.547} // K & L
    };
    private const double L23 = 0.905129, L34 = 1.308227;
    private readonly double[] algaeHeights = new double[] {L23, L34};
    private readonly double[,] blueAlgae = new double[,] {
        {3.803, 4.025, L34}, // A & B
        {4.144, 3.42, L23}, // C & D
        {4.836, 3.42, L34}, // E & F
        {5.177, 4.025, L23}, // G & H
        {4.836, 4.61, L34}, // I & J
        {4.144, 4.61, L23} // K & L
    };
    private readonly double[,] redAlgae = new double[,] {
        {13.737, 4.025, L34}, // A & B
        {13.396, 4.61, L23}, // C & D
        {12.704, 4.61, L34}, // E & F
        {12.363, 4.025, L23}, // G & H
        {12.704, 3.42, L34}, // I & J
        {13.396, 3.42, L23} // K & L
    };

    private double[,] coralXYPositions = new double[0, 0], algaePositions = new double[0, 0];

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        // Sets up network tables
        // TODO: Change serverBaseAddress and potentially add serverPort, onOpen, and onNewTopicData
        _ntClient = new Nt4Client("Reefbreaker", "127.0.0.1", onNewTopicData: OnNewTopicData);
        _ntClient.Connect();

        // Subscribes to relevant topics
        _cameraOffsetsSubscriptionId = _ntClient.Subscribe("cameraOffsets");
        _coralCameraAnglesSubscriptionId = _ntClient.Subscribe("coralCameraAngles");
        _algaeCameraAnglesSubscriptionId = _ntClient.Subscribe("coralCameraAngles");
        _robotPositionSubscriptionId = _ntClient.Subscribe("robotPosition");
        _allianceSubscriptionId = _ntClient.Subscribe("alliance");
        _cameraParametersSubscriptionId = _ntClient.Subscribe("cameraParameters");
    }

    // Update is called once per frame
    void Update()
    {
        coralXYPositions = alliance == "Blue" ? blueCoral : redCoral;
        algaePositions = alliance == "Blue" ? blueAlgae : redAlgae;
        updateCameraPosition();
        coralRayConstants = RayConstantsGenerator(coralCameraAngles);
        algaeRayConstants = RayConstantsGenerator(algaeCameraAngles);
        coralRawIntersections = RawIntersectionFinder(coralRayConstants);
        algaeRawIntersections = RawIntersectionFinder(algaeRayConstants);
        coralDistanceEstimates = getDistanceEstimates(coralCameraAngles);
        algaeDistanceEstimates = getDistanceEstimates(algaeCameraAngles);
        coralValidIntersections = declareValidIntersections(coralRawIntersections, coralDistanceEstimates);
        algaeValidIntersections = declareValidIntersections(algaeRawIntersections, algaeDistanceEstimates);
        coralCorrectIntersections = findCorrectIntersectionsAndUpdateCounters(coralValidIntersections, hasSeenCoral, "coral");
        algaeCorrectIntersections = findCorrectIntersectionsAndUpdateCounters(algaeValidIntersections, hasSeenAlgae, "algae");
        coralVisibilities = checkCoralVisibilies();
        algaeVisibilities = checkAlgaeVisibilies();
        updateEmptyCounters(emptyCoralCounters, hasSeenCoral, coralVisibilities);
        updateEmptyCounters(emptyAlgaeCounters, hasSeenAlgae, algaeVisibilities);
        checkEmptyCoralCounters();
        updateGamepiecePresence(coralPresence, coralCounters);
        updateGamepiecePresence(algaePresence, algaeCounters);
        publishValues();
        resetValues();
    }

    // Method run each time new values are posted to network tables
    // TODO: Ensure robot code publishes jagged arrays for cameraAngles or change method
    private void OnNewTopicData(Nt4Topic topic, long timestamp, object value)
    {
        if (topic.Name == "robotPosition" && value is double[] robotPose && robotPose.Length == 4)
        {
            robotX = robotPose[0];
            robotY = robotPose[1];
            robotZ = robotPose[2];
            robotRotation = robotPose[3] + 180; // TODO: Check assumption that robot rotation is from -180 to 180
        }
        else if (topic.Name == "coralCameraAngles" && value is double[][][] coralCameraAngles) // cameraAngles[camera][angles][xAngle, yAngle, width, height]
        {
            this.coralCameraAngles = coralCameraAngles;
        }
        else if (topic.Name == "algaeCameraAngles" && value is double[][][] algaeCameraAngles) // cameraAngles[camera][angles][xAngle, yAngle, width, height]
        {
            this.algaeCameraAngles = algaeCameraAngles;
        }
        else if (topic.Name == "alliance" && value is String alliance) // "Blue" or "Red"
        {
            this.alliance = alliance;
        }
        else if (topic.Name == "cameraOffsets" && value is double[][] cameraOffsets) // cameraOffset[xOffset, yOffset, zOffset, pitch]
        {
            this.cameraOffsets = cameraOffsets;
        }
        else if (topic.Name == "cameraParameters" && value is int[][] cameraParameters) // cameraParameters[camera][xResolution, yResolution, xFov, yFov]
        {
            this.cameraParameters = cameraParameters; 
        }
    }

    // Likely will not need this function as camera angles will probably be returned
    private double[][][] correctCameraAngles(double[][][] cameraAngles)
    {
        double[][][] correctAngles = copy3DArray(cameraAngles); // [camera][angles][xPixel -> xAngle, yPixel -> yAngle, width, height]
        int xPixel, yPixel, xCenter, yCenter;
        double xFocal, yFocal, xAngle, yAngle;
        for (int camera = 0; camera < correctAngles.Length; camera++)
        {
            for (int angle = 0; angle < correctAngles[camera].Length; angle++)
            {
                double[] angleParameters = correctAngles[camera][angle]; // [xPixel, yPixel, width, height]
                int[] cameraParameters = this.cameraParameters[camera]; // [xResolution, yResolution, xFov, yFov]
                xCenter = cameraParameters[0] / 2;
                yCenter = cameraParameters[1] / 2;
                xPixel = (int)angleParameters[0] - xCenter;
                yPixel = (int)angleParameters[1] - yCenter;
                xFocal = angleParameters[2] / (2 * Math.Tan(cameraParameters[2] / 2.0));
                yFocal = angleParameters[3] / (2 * Math.Tan(cameraParameters[3] / 2.0));
                xAngle = Math.Atan(xPixel / xFocal);
                yAngle = Math.Atan(yPixel / yFocal);
                correctAngles[camera][angle][0] = xAngle;
                correctAngles[camera][angle][1] = yAngle;
            }
        }
        return correctAngles;
    }

    // returns cameraPositions[camera][x, y, z, pitch]
    private void updateCameraPosition()
    {
        cameraPositions = copy2DArraySize(cameraOffsets);
        for (int camera = 0; camera < cameraOffsets.Length; camera++)
        {
            cameraPositions[camera][0] = robotX + cameraOffsets[camera][0] * Math.Cos(robotRotation * (Math.PI / 180)); // X
            cameraPositions[camera][1] = robotY + cameraOffsets[camera][1] * Math.Sin(robotRotation * (Math.PI / 180)); // Y
            cameraPositions[camera][2] = cameraOffsets[camera][2]; // Z
            cameraPositions[camera][3] = cameraOffsets[camera][3]; // Pitch
            cameraPositions[camera][4] = cameraOffsets[camera][4]; // Yaw
        }
    }

    // TODO: Check assumption that cameraX is left-/right+ while cameraY is up+/down- in frame
    // Returns rayConstants[camera][angle][constant]
    private double[][][] RayConstantsGenerator(double[][][] cameraAngles)
    {
        double[][][] rayConstants = copyOuter2Dof3DArrayWithFinalColumn(cameraAngles, 6); // 0=ax, 1=mx, 2=ay, 3=my, 4=az, 5=mz
        for (int camera = 0; camera < cameraAngles.Length; camera++)
        {
            for (int angle = 0; angle < cameraAngles[camera].Length; angle++)
            {
                rayConstants[camera][angle][0] = cameraPositions[camera][0];
                rayConstants[camera][angle][1] = Math.Cos((cameraAngles[camera][angle][1] + cameraPositions[camera][3]) * (Math.PI / 180)); // Cos(yAngle + pitch)
                rayConstants[camera][angle][2] = cameraPositions[camera][1];
                rayConstants[camera][angle][3] = Math.Sin(-cameraAngles[camera][angle][0] * (Math.PI / 180)); // Sin(-xAngle)
                rayConstants[camera][angle][4] = cameraPositions[camera][2];
                rayConstants[camera][angle][5] = Math.Sin((cameraAngles[camera][angle][1] + cameraPositions[camera][3]) * (Math.PI / 180)); // Sin(yAngle + pitch)
            }
        }
        return rayConstants;
    }

    // Returns intersections[intersection set][intersection1, intersection2][x, y, z]
    private double[][][] RawIntersectionFinder(double[][][] rayConstants)
    {
        double[][][] rawIntersections = createNew3DArray(getNumOf1DArraysOf3DArray(rayConstants), 2, 3);
        double a, b, c, reefX, parameter1, parameter2;
        double[] constants;
        int count = 0;
        for (int camera = 0; camera < rayConstants.Length; camera++)
        {
            for (int angle = 0; angle < rayConstants[camera].Length; angle++)
            {
                // Reference: 0=ax, 1=mx, 2=ay, 3=my, 4=az, 5=mz
                constants = rayConstants[camera][angle];
                reefX = alliance == "Blue" ? blueReefX : reedReefX;
                a = Math.Pow(reefX, 2) + Math.Pow(reefY, 2); // (mx)^2 + (my)^2
                b = 2 * (constants[1] * (constants[0] - reefX) + constants[3] * (constants[2] - reefY)); // 2((mx)(ax - reefX) + (my)(ay - reefY))
                c = Math.Pow(constants[0] - reefX, 2) + Math.Pow(constants[2] - reefY, 2) - Math.Pow(radius, 2); // (ax - reefX)^2 + (ay - reefY)^2 - r^2
                parameter1 = (-b + Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);
                parameter2 = (-b - Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);
                rawIntersections[count][0] = new double[] { constants[0] + (constants[1] * parameter1), constants[2] + (constants[3] * parameter1), constants[4] + (constants[5] * parameter1) };
                rawIntersections[count][1] = new double[] { constants[0] + (constants[1] * parameter2), constants[2] + (constants[3] * parameter2), constants[4] + (constants[5] * parameter2) };
                count++;
            }
        }
        return rawIntersections;
    }

    private double getEstimatedDistance(double width, double height)
    {
        // TODO: Create a function based on collected data to estimate the distance
        return 0.0;
    }

    // Returns estimates[estimate]
    private double[] getDistanceEstimates(double[][][] cameraAngles)
    {
        double[] distanceEstimates = new double[getNumOf1DArraysOf3DArray(cameraAngles)];
        int count = 0;
        for (int camera = 0; camera < cameraAngles.Length; camera++)
        {
            for (int angle = 0; angle < cameraAngles[camera].Length; angle++)
            {
                distanceEstimates[count] = getEstimatedDistance(cameraAngles[camera][angle][2], cameraAngles[camera][angle][3]);
                count++;
            }
        }
        return distanceEstimates;
    }

    // returns intersections[intersection][x, y, z]
    private double[][] declareValidIntersections(double[][][] rawIntersections, double[] distanceEstimates)
    {
        double[][] validIntersections = createNew2DArray(distanceEstimates.Length, 3);
        double[] robotPosition = new double[] { robotX, robotY, robotZ };
        double distance1, distance2;
        for (int intersectionNum = 0; intersectionNum < rawIntersections.Length; intersectionNum++)
        {
            distance1 = get3DDistance(rawIntersections[intersectionNum][0], robotPosition);
            distance2 = get3DDistance(rawIntersections[intersectionNum][1], robotPosition);
            validIntersections[intersectionNum] = distance1 - distanceEstimates[intersectionNum] < distance2 - distanceEstimates[intersectionNum] ? rawIntersections[intersectionNum][0] : rawIntersections[intersectionNum][1];
        }
        return validIntersections;
    }

    // returns points[point][x, y, z]
    private double[][] findCorrectIntersectionsAndUpdateCounters(double[][] validIntersections, bool[] hasSeen, String gamepiece) {
        double[][] correctIntersections = copy2DArraySize(validIntersections);
        double[,] correctXYs = coralXYPositions;
        double[] correctZs = new double[0];
        int[] counters = new int[0];
        if (gamepiece == "coral") {
            correctZs = coralHeights;
            counters = coralCounters;
        } else if (gamepiece == "algae") {
            correctZs = algaeHeights;
            counters = algaeCounters;
        }
        for (int i = 0; i < correctIntersections.Length; i++) {
            correctIntersections[i] = findCorrectIntersectionAndUpdateCounter(validIntersections[i], correctXYs, correctZs, counters, hasSeen, gamepiece);
        }
        return correctIntersections;
    }

    // returns points[x, y, z]
    private double[] findCorrectIntersectionAndUpdateCounter(double[] validIntersection, double[,] correctXYs, double[] correctZs, int[] counters, bool[] hasSeen, String gamepiece) {
        double correctX = 0, correctY = 0, correctZ = 0, smallestXYDistance = 10, smallestZDistance = 10, xyDistance, zDistance;
        int correctXYIndex = 0, correctZIndex = 0, correctIndex = 0;
        for (int i = 0; i < correctXYs.GetLength(0); i++) {
            double[] XY = new double[] {correctXYs[i, 0], correctXYs[i, 1]};
            xyDistance = get2DDistance(validIntersection, XY);
            if (xyDistance < smallestXYDistance) {
                correctX = XY[0];
                correctY = XY[1];
                smallestXYDistance = xyDistance;
                correctXYIndex = i;
            }
        }
        for (int i = 0; i < correctZs.Length; i++) {
            zDistance = Math.Abs(validIntersection[2] - correctZs[i]);
            if (zDistance < smallestZDistance) {
                correctZ = correctZs[i];
                smallestZDistance = zDistance;
                correctZIndex = i;
            }
        }
        correctIndex = gamepiece == "coral" ?  correctXYIndex * correctZs.Length + correctZIndex : correctXYIndex;
        counters[correctIndex]++;
        hasSeen[correctIndex] = true;
        return new double[] {correctX, correctY, correctZ};
    }

    private void checkEmptyCoralCounters() {
        for (int i = 0; i < emptyCoralCounters.Length; i++) {
            if (emptyCoralCounters[i] > 0) {
                coralCounters[i] = 0;
            }
        }
    }

    private void updateGamepiecePresence(bool[] pieceStates, int[] counters) {
        for (int i = 0; i < counters.Length; i++) {
            if (counters[i] > framePresenceThreshold) {
                pieceStates[i] = true;
            } else if (counters[i] == 0) {
                pieceStates[i] = false;
            }
        }
    }

    private void updateEmptyCountersBasedOnVisibility(double[] emptyCounters, double[,] correctXYs, double[] correctZs) {
        for (int i = 0; i < correctXYs.Length; i++) {
            for (int j = 0; j < correctZs.Length; j++) {

            }
        }
    }

    // returns visible[coral]
    private bool[] checkCoralVisibilies() {
        bool[] coralVisibilities = new bool[coralCounters.Length];
        for (int i = 0; i < cameraPositions.Length; i++) {
            for (int j = 0; j < coralXYPositions.GetLength(0); j++) {
                for (int k = 0; k < coralHeights.Length; k++) {
                    double[] cameraPosition = cameraPositions[i];
                    if (IsVisible(new Vector3((float) cameraPosition[0], (float) cameraPosition[1], (float) cameraPosition[2]), (float) cameraPosition[4], (float) cameraPosition[3], (float) cameraParameters[i][2], (float) cameraParameters[i][3], new Vector3((float) coralXYPositions[j, 0], (float) coralXYPositions[j, 1], (float) coralHeights[k]))) {
                        coralVisibilities[j * coralHeights.Length + k] = true;
                    }
                }
            }
        }
        return coralVisibilities;
    }

    private bool[] checkAlgaeVisibilies() {
        bool[] algaeVisibilities = new bool[algaeCounters.Length];
        for (int i = 0; i < cameraPositions.Length; i++) {
            for (int j = 0; j < algaePositions.GetLength(0); j++) {
                double[] cameraPosition = cameraPositions[i];
                    if (IsVisible(new Vector3((float) cameraPosition[0], (float) cameraPosition[1], (float) cameraPosition[2]), (float) cameraPosition[4], (float) cameraPosition[3], (float) cameraParameters[i][2], (float) cameraParameters[i][3], new Vector3((float) algaePositions[j, 0], (float) algaePositions[j, 1], (float) algaePositions[j, 2]))) {
                        algaeVisibilities[j] = true;
                    }
            }
        }
        return algaeVisibilities;
    }

    private bool IsVisible(Vector3 cameraPos, float yaw, float pitch, float fovX, float fovY, Vector3 objectPos)
    {
        // Convert degrees to radians
        yaw = Mathf.Deg2Rad * yaw;
        pitch = Mathf.Deg2Rad * pitch;
        float halfFovX = fovX * Mathf.Deg2Rad / 2;
        float halfFovY = fovY * Mathf.Deg2Rad / 2;

        // Compute camera coordinate system (right, up, forward vectors)
        Vector3 forward = new Vector3(
            Mathf.Cos(pitch) * Mathf.Cos(yaw),
            Mathf.Sin(pitch),
            Mathf.Cos(pitch) * Mathf.Sin(yaw)
        );

        Vector3 right = new Vector3(-Mathf.Sin(yaw), 0, Mathf.Cos(yaw));  // Right direction
        Vector3 up = Vector3.Cross(forward, right);  // Up direction using cross product

        // Construct the transformation matrix
        Matrix4x4 cameraMatrix = new Matrix4x4();
        cameraMatrix.SetColumn(0, new Vector4(right.x, right.y, right.z, 0));
        cameraMatrix.SetColumn(1, new Vector4(up.x, up.y, up.z, 0));
        cameraMatrix.SetColumn(2, new Vector4(forward.x, forward.y, forward.z, 0));
        cameraMatrix.SetColumn(3, new Vector4(cameraPos.x, cameraPos.y, cameraPos.z, 1));

        // Convert object position to camera space
        Vector3 relativePos = objectPos - cameraPos;
        Vector3 objectInCameraSpace = cameraMatrix.inverse.MultiplyPoint3x4(relativePos);

        // Get the transformed coordinates
        float Vx = objectInCameraSpace.x;
        float Vy = objectInCameraSpace.y;
        float Vz = objectInCameraSpace.z;

        // Ensure the object is in front of the camera
        if (Vz <= 0) return false;

        // Compute angles
        float thetaX = Mathf.Atan2(Vx, Vz);
        float thetaY = Mathf.Atan2(Vy, Vz);

        // Check if within FOV
        return Mathf.Abs(thetaX) <= halfFovX && Mathf.Abs(thetaY) <= halfFovY;
    }

    private void updateEmptyCounters(int[] emptyCounters, bool[] hasSeen, bool[] isVisible) {
        for (int i = 0; i < emptyCounters.Length; i++) {
            if (isVisible[i] && !hasSeen[i]) {
                emptyCounters[i]++;
            } else if (hasSeen[i]) {
                emptyCounters[i] = 0;
            }
        }
    }

    private void publishValues() {
        _ntClient.PublishValue("Coral Presence", coralPresence);
        _ntClient.PublishValue("Algae Presence", algaePresence);
    }

    private void resetValues() {
        hasSeenCoral = new bool[12];
        hasSeenAlgae = new bool[12];
        coralVisibilities = new bool[0];
        algaeVisibilities = new bool[0];
        coralDistanceEstimates = new double[0];
        algaeDistanceEstimates = new double[0];
        coralValidIntersections = new double[0][];
        algaeValidIntersections = new double[0][];
        coralCorrectIntersections = new double[0][];
        algaeCorrectIntersections = new double[0][];
        coralCameraAngles = new double[0][][];
        algaeCameraAngles = new double[0][][];
        coralRayConstants = new double[0][][];
        algaeRayConstants = new double[0][][];
        coralRawIntersections = new double[0][][];
        algaeRawIntersections = new double[0][][];
    }

    private double[][] copy2DArraySize(double[][] original)
    {
        double[][] newArray = new double[original.Length][];
        for (int i = 0; i < original.Length; i++)
        {
            newArray[i] = new double[original[i].Length];
        }
        return newArray;
    }

    private double[][][] copyOuter2Dof3DArrayWithFinalColumn(double[][][] original, int lastDimension)
    {
        double[][][] newArray = new double[original.Length][][];
        for (int i = 0; i < original.Length; i++)
        {
            newArray[i] = new double[original[i].Length][];
            for (int j = 0; j < original[i].Length; j++)
            {
                newArray[i][j] = new double[lastDimension];
            }
        }
        return newArray;
    }

    private double[][][] createNew3DArray(int D1, int D2, int D3)
    {
        double[][][] array = new double[D1][][];
        for (int i = 0; i < D1; i++)
        {
            array[i] = new double[D2][];
            for (int j = 0; j < D2; j++)
            {
                array[i][j] = new double[D3];
            }
        }
        return array;
    }

    private double[][] createNew2DArray(int D1, int D2)
    {
        double[][] newArray = new double[D1][];
        for (int i = 0; i < D1; i++)
        {
            newArray[i] = new double[D2];
        }
        return newArray;
    }

    private int getNumOf1DArraysOf3DArray(double[][][] array)
    {
        int numArrays = 0;
        for (int i = 0; i < array.Length; i++)
        {
            numArrays += array[i].Length;
        }
        return numArrays;
    }

    private double[][][] copy3DArray(double[][][] original)
    {
        double[][][] newArray = new double[original.Length][][];
        for (int i = 0; i < original.Length; i++)
        {
            newArray[i] = new double[original[i].Length][];
            for (int j = 0; j < original[i].Length; j++)
            {
                newArray[i][j] = new double[original[i][j].Length];
                Array.Copy(original[i][j], newArray[i][j], original[i][j].Length);
            }
        }
        return newArray;
    }

    private double get3DDistance(double[] point1, double[] point2)
    {
        return Math.Sqrt(Math.Pow(point2[0] - point1[0], 2) + Math.Pow(point2[1] - point1[1], 2) + Math.Pow(point2[2] - point1[2], 2));
    }

    private double get2DDistance(double[] point1, double[] point2) {
        return Math.Sqrt(Math.Pow(point2[0] - point1[0], 2) + Math.Pow(point2[1] - point1[1], 2));
    }
}