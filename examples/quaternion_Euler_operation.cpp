#include <gsCore/gsQuaternion.h>
#include <iostream>
#include <gismo.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef GISMO_DEG2RAD
#define GISMO_DEG2RAD (M_PI / 180.0)
#endif

#ifndef GISMO_RAD2DEG
#define GISMO_RAD2DEG (180.0 / M_PI)
#endif

using namespace gismo;


int main(int argc, char* argv[])
{
    // Create a quaternion using the default constructor
    gsQuaternion<real_t> q1;
    gsVector<real_t> euler(3);
    gsInfo << "Default quaternion: [" << q1.e0() << ", " << q1.e1() << ", " << q1.e2() << ", " << q1.e3() << "]" << "\n";

    double alpha1 = 10.0;
    double beta1  = 11.0;
    double gamma1 = 12.0;

    gsInfo << "Rotation about frame axes:" << "\n";

    // Test rotation about X-axis
    gsInfo << "Rotation about x-axis: " << alpha1 << " degrees" << "\n";
    q1.SetFromAngleX(GISMO_DEG2RAD * alpha1);
    euler = GISMO_RAD2DEG * q1.GetCardanXyz();
    gsInfo << "Euler angles: [" << euler << "]" << "\n";

    // Test rotation about Y-axis
    gsInfo << "Rotation about y-axis: " << beta1 << " degrees" << "\n";
    q1.SetFromAngleY(GISMO_DEG2RAD * beta1);
    euler = GISMO_RAD2DEG * q1.GetCardanXyz();
    gsInfo << "Euler angles: [" << euler << "]" << "\n";

    // Test rotation about Z-axis
    gsInfo << "Rotation about z-axis: " << gamma1 << " degrees" << "\n";
    q1.SetFromAngleZ(GISMO_DEG2RAD * gamma1);
    euler = GISMO_RAD2DEG * q1.GetCardanXyz();
    gsInfo << "Euler angles: [" << euler << "]" << "\n";

    //-------------------------------Generate quaternion from Euler angles--------------------------------
    gsInfo << "Quaternion from Euler angles:" << "\n";

    euler << alpha1, beta1, gamma1;
    gsInfo << "Euler angles: [" << euler << "]" << "\n";
    // Convert Euler angles (in degrees) to radians and create quaternion
    q1.SetCardanXyz(GISMO_DEG2RAD * euler);
    gsInfo << "Quaternion: [" << q1.e0() << ", " << q1.e1() << ", " << q1.e2() << ", " << q1.e3() << "]" << "\n";
    euler = GISMO_RAD2DEG * q1.GetCardanXyz();
    gsInfo << "Euler angles (recovered): [" << euler << "]" << "\n";

    // Test with another set of Euler angles
    real_t alpha2 = -17.3;
    real_t beta2  = -41.0;
    real_t gamma2 = 0.7;
    gsVector<real_t> euler2(3);
    gsQuaternion<real_t> q2;
    euler2 << alpha2, beta2, gamma2;
    gsInfo << "Euler angles: [" << euler2 << "]" << "\n";
    q2.SetCardanXyz(GISMO_DEG2RAD * euler2);
    gsInfo << "Quaternion: [" << q2.e0() << ", " << q2.e1() << ", " << q2.e2() << ", " << q2.e3() << "]" << "\n";
    euler2 = GISMO_RAD2DEG * q2.GetCardanXyz();
    gsInfo << "Euler angles (recovered): [" << euler2 << "]" << "\n";

    //-------------------------------Rotation vector and unit quaternion test--------------------------------
    gsQuaternion<real_t> q_unit;
    q_unit.SetUnit();
    std::cout << "Unit quaternion (constructed from scalar): ["
              << q_unit.e0() << ", " << q_unit.e1() << ", "
              << q_unit.e2() << ", " << q_unit.e3() << "]" << std::endl;

    gsInfo << "q_unit is unit quaternion: " << q_unit.IsIdentity() << "\n";

    // Test rotation of a vector using q_unit (which represents identity, so no rotation)
    gsVector<real_t> vec(3);
    vec << 1, 2, 3;
    gsVector<real_t> vec_rot = q_unit.Rotate(vec);
    gsVector<real_t> vec_back = q_unit.RotateBack(vec_rot);
    std::cout << "Original vector: " << vec << std::endl;
    std::cout << "Rotated vector: " << vec_rot << std::endl;
    std::cout << "Recovered vector (after inverse rotation): " << vec_back << std::endl;

    // -------------------- Test derivative calculations --------------------
    // Suppose an angular speed (absolute coordinates)
    gsVector<real_t> ang_speed(3);
    ang_speed << 0.1, 0.2, 0.3;
    gsQuaternion<real_t> q_dt;
    q_dt.SetDerivativeAbs(ang_speed, q_unit);
    std::cout << "First derivative of q (absolute): ["
              << q_dt.e0() << ", " << q_dt.e1() << ", " << q_dt.e2() << ", " << q_dt.e3() << "]" << std::endl;

    // Suppose an angular acceleration (absolute coordinates)
    gsVector<real_t> ang_acc(3);
    ang_acc << 0.01, 0.02, 0.03;
    gsQuaternion<real_t> q_dtdt;
    q_dtdt.SetDt2FromAngAccAbs(ang_acc, q_unit, q_dt);
    std::cout << "Second derivative of q (absolute): ["
              << q_dtdt.e0() << ", " << q_dtdt.e1() << ", " << q_dtdt.e2() << ", " << q_dtdt.e3() << "]" << std::endl;

    // Test GetAngularSpeedAbs
    gsVector<real_t> computed_speed(3);
    q_dt.GetAngularSpeedAbs(computed_speed, q_unit);
    std::cout << "Angular speed (absolute) computed from q_dt: " << computed_speed << std::endl;

    // -------------------- Test relative derivative calculations --------------------
    gsVector<real_t> ang_speed_rel(3);
    ang_speed_rel << 0.05, 0.06, 0.07;
    gsQuaternion<real_t> q_dt_rel;
    q_dt_rel.SetDerivativeRel(ang_speed_rel, q_unit);
    std::cout << "First derivative of q (relative): ["
              << q_dt_rel.e0() << ", " << q_dt_rel.e1() << ", " << q_dt_rel.e2() << ", " << q_dt_rel.e3() << "]" << std::endl;
    gsVector<real_t> ang_acc_rel(3);
    ang_acc_rel << 0.005, 0.006, 0.007;
    gsQuaternion<real_t> q_dtdt_rel;
    q_dtdt_rel.SetDt2FromAngAccRel(ang_acc_rel, q_unit, q_dt_rel);
    std::cout << "Second derivative of q (relative): ["
              << q_dtdt_rel.e0() << ", " << q_dtdt_rel.e1() << ", " << q_dtdt_rel.e2() << ", " << q_dtdt_rel.e3() << "]" << std::endl;
    
    // Test GetAngularSpeedRel
    gsVector<real_t> computed_speed_rel(3);
    q_dt_rel.GetAngularSpeedRel(computed_speed_rel, q_unit);
    std::cout << "Angular speed (relative) computed from q_dt_rel: " << computed_speed_rel << std::endl;

    return 0;
}