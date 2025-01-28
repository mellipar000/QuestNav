using System;
using UnityEngine;

public class RayTracing : MonoBehaviour
{
    // Angles to game pieces relative to camera
    private double xAngle, yAngle, zAngle;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public double[] getCameraAngles() {
        xAngle = 0.0;
        yAngle = 0.0;
        zAngle = 0.0;
        return new double[] {xAngle, yAngle, zAngle};
    }
}
