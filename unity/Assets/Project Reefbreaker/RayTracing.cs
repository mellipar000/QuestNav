using System;
using NetworkTables;
using UnityEngine;

public class RayTracing : MonoBehaviour
{
    private Nt4Client _ntClient;
    private double robotX = 0.0, robotY = 0.0, robotZ = 0.0, robotRotation = 0.0;
    private double[,] cameraOffsets = new double[0, 0], cameraPositions = new double[0, 0], parameters = new double[0, 0];
    private double[][][] cameraAngles = new double[0][][], rayConstants = new double[0][][];
    private String alliance = "Blue";
    private int _cameraOffsetSubscriptionId, _cameraAnglesSubscriptionId, _robotPositionSubscriptionId, _allianceSubscriptionId;

    // Field Constants
    private const double blueReefX = 4.48945, reedReefX = 13.065, reefY = 4.0259, radius = 0.658;
    private const double L1 = 0.4572, L2 = 0.70358, L3 = 1.106678, L4 = 1.652778;
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

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        // Sets up network tables
        // TODO: Change serverBaseAddress and potentially add serverPort, onOpen, and onNewTopicData
        _ntClient = new Nt4Client("Reefbreaker", "127.0.0.1", onNewTopicData: OnNewTopicData);
        _ntClient.Connect();

        // Subscribes to relevant topics
        _cameraOffsetSubscriptionId = _ntClient.Subscribe("cameraOffsets");
        _cameraAnglesSubscriptionId = _ntClient.Subscribe("cameraAngles");
        _robotPositionSubscriptionId = _ntClient.Subscribe("robotPosition");
        _allianceSubscriptionId = _ntClient.Subscribe("alliance");
    }

    // Update is called once per frame
    void Update()
    {
        updateCameraPosition();
        RayConstantsGenerator();
        ParameterFinder();
    }

    // Method run each time new values are posted to network tables
    // TODO: Ensure robot code publishes jagged arrays for cameraAngles or change method
    private void OnNewTopicData(Nt4Topic topic, long timestamp, object value)
    {
        if (topic.Name == "robotPosition" && value is double[] robotPose && robotPose.Length == 3)
        {
            robotX = robotPose[0];
            robotY = robotPose[1];
            robotZ = robotPose[2];
        }
        else if (topic.Name == "cameraAngles" && value is double[][][] cameraAngles)
        {
            this.cameraAngles = cameraAngles;
        }
        else if (topic.Name == "alliance" && value is String alliance)
        {
            this.alliance = alliance;
        }
        else if (topic.Name == "cameraOffsets" && value is double[,] cameraOffsets)
        {
            this.cameraOffsets = cameraOffsets;
        }
    }

    private void updateCameraPosition()
    {
        cameraPositions = new double[cameraOffsets.GetLength(0), cameraOffsets.GetLength(1)];
        // TODO: Ensure robot rotation goes from 0 to 360, not -180 to 180
        for (int camera = 0; camera < cameraOffsets.GetLength(0); camera++)
        {
            cameraPositions[camera, 0] = robotX + cameraOffsets[camera, 0] * Math.Cos(robotRotation * (Math.PI / 180)); // X
            cameraPositions[camera, 1] = robotY + cameraOffsets[camera, 1] * Math.Sin(robotRotation * (Math.PI / 180)); // Y
            cameraPositions[camera, 2] = cameraOffsets[camera, 2]; // Z
            cameraPositions[camera, 3] = cameraOffsets[camera, 3]; // Pitch
        }
    }

    // TODO: Check assumption that cameraX is left-/right+ while cameraY is up+/down- in frame
    private void RayConstantsGenerator()
    {
        rayConstants = Create3DJaggedArrayWithSameBase(cameraAngles, 6); // 0=ax, 1=mx, 2=ay, 3=my, 4=az, 5=mz
        for (int camera = 0; camera < cameraAngles.Length; camera++)
        {
            for (int angle = 0; angle < cameraAngles[camera].Length; angle++)
            {
                rayConstants[camera][angle][0] = cameraPositions[camera, 0];
                rayConstants[camera][angle][1] = Math.Cos((cameraAngles[camera][angle][1] + cameraPositions[camera, 3]) * (Math.PI / 180)); // Cos(yAngle + pitch)
                rayConstants[camera][angle][2] = cameraPositions[camera, 1];
                rayConstants[camera][angle][3] = Math.Sin(-cameraAngles[camera][angle][0] * (Math.PI / 180)); // Sin(-xAngle)
                rayConstants[camera][angle][4] = cameraPositions[camera, 2];
                rayConstants[camera][angle][5] = Math.Sin((cameraAngles[camera][angle][1] + cameraPositions[camera, 3]) * (Math.PI / 180)); // Sin(yAngle + pitch)
            }
        }
    }

    private void ParameterFinder()
    {
        parameters = new double[GetTotalNumArray3dJagged(rayConstants), 2];
        double a, b, c, reefX;
        int count = 0;
        double[] constants;
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
                parameters[count, 0] = (-b + Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);
                parameters[count, 1] = (-b - Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);
                count++;
            }
        }
    }

    private static double[][][] Create3DJaggedArrayWithSameBase(double[][][] original, int finalColumnNumber)
    {
        double[][][] newArray = new double[original.Length][][];
        for (int i = 0; i < original.Length; i++)
        {
            newArray[i] = new double[original[i].Length][];
            for (int j = 0; j < newArray[i].Length; j++)
            {
                newArray[i][j] = new double[finalColumnNumber];
            }
        }
        return newArray;
    }

    public int GetTotalNumArray3dJagged(double[][][] jaggedArray)
    {
        int count = 0;
        for (int i = 0; i < jaggedArray.Length; i++)
        {
            count += jaggedArray[i].Length;
        }
        return count;
    }
}
