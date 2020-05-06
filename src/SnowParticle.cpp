#include "SnowParticle.hpp"

SnowParticleMaterial::SnowParticleMaterial()
{
    mu = youngsModule / (2. + 2. * PoissonsRatio);
    lambda = youngsModule * PoissonsRatio /
             ((1. + PoissonsRatio) * (1. - 2. * PoissonsRatio));
}

SnowParticleMaterial::~SnowParticleMaterial()
{
}

SnowParticle::SnowParticle() : m(nullptr)
{
    // TODO proper constructor
}

SnowParticle::SnowParticle(const Vector3f& pos, SnowParticleMaterial* material)
    : position(pos), m(material)
{
    // TODO proper constructor
}

SnowParticle::~SnowParticle()
{
}

void SnowParticle::updatePos()
{
    // Simple euler integration
    position += deltaT * velocity;
}
void SnowParticle::updatePureElasticGradient()
{
    // calculate the pure elastic deformation after this time step
    Matrix3f pureEDefGradChangeRate =
        velocityGradient * deltaT + Matrix3f::Identity();
    deformationGradientElastic =
        pureEDefGradChangeRate * deformationGradientElastic;
}
void SnowParticle::updateCombinedPElasticGradient()
{
    // TODO very careful about row/col
    // TODO make sure about the row/or col major is consistent with the paper
    // and with eigen

    // consider both elastic and plastic then
    Matrix3f f_all = deformationGradientElastic * deformationGradientPlastic;
    // We compute the SVD decomposition
    // The singular values (basically a scale transform) tell us if
    // the particle has exceeded critical stretch/compression

    // TODO find out how eigen do SVD decompesation
    // TODO
    // deformationGradientElastic.svd(&SVDW, &SVDE, &SVDV);
    Matrix3f SVDV_trans = SVDV.transpose();
    // Clamp singular values to within elastic region
    for (int i = 0; i < 3; i++)
    {
        // if (SVDE[i] < m->criticalStress)
        //     SVDE[i] = m->criticalStress;
        // else if (SVDE[i] > m->criticalStretch)
        //     SVDE[i] = m->criticalStretch;

        SVDE[i] = clamp(m->criticalStress, m->criticalStretch, SVDE[i]);
    }

    // TODO fig out if 3d needs polar-phi, how to compute
    // TODO
    // #if ENABLE_IMPLICIT
    // Compute polar decomposition, from clamped SVD
    // polar_r.setData(SVDW * SVDV_trans);
    // polar_s.setData(SVDV);
    // polar_s.diag_product(SVDE);
    // polar_s.setData(polar_s * SVDV_trans);
    // #endif

    // Recompute elastic and plastic gradient
    // We're basically just putting the SVD back together again

    Matrix3f vDivideSVED(SVDV), wTimeSVDE(SVDW);
    vDivideSVED.row(0) /= SVDE[0];
    vDivideSVED.row(1) /= SVDE[1];
    vDivideSVED.row(2) /= SVDE[2];
    wTimeSVDE.row(0) *= SVDE[0];
    wTimeSVDE.row(1) *= SVDE[1];
    wTimeSVDE.row(2) *= SVDE[2];
    // v_cpy.diag_product_inv(SVDE);
    // w_cpy.diag_product(SVDE);
    deformationGradientPlastic = vDivideSVED * SVDW.transpose() * f_all;
    deformationGradientElastic = wTimeSVDE * SVDV.transpose();
}
const Matrix3f SnowParticle::energyDerivative()
{
    // Adjust lame parameters to account for m->hardening
    float harden =
        exp(m->hardening * (1. - deformationGradientPlastic.determinant()));
    float Je = SVDE.x() * SVDE.y() * SVDE.z();
    // This is the co-rotational term
    Matrix3f temp = 2. * m->mu *
                    (deformationGradientElastic - SVDW * SVDV.transpose()) *
                    deformationGradientElastic.transpose();
    // Add in the primary contour term
    temp += m->lambda * Je * (Je - 1.) * Matrix3f::Identity();
    // temp.diag_sum(m->lambda * Je * (Je - 1));
    // Add m->hardening and volume
    return volume * harden * temp;
}
// #if ENABLE_IMPLICITF
const Vector3f SnowParticle::deltaForce(const Vector2f& u,
                                        const Vector2f& weight_grad)
{
    return Vector3f(1., 1., 1.);
    // TODO

    // // For detailed explanation, check out the implicit math pdf for details
    // // Before we do the force calculation, we need deltaF, deltaR, and
    // // delta(JF^-T)

    // // Finds delta(Fe), where Fe is the elastic deformation gradient
    // // Probably can optimize this expression with parentheses...
    // Matrix2f del_elastic =
    //     deltaT * u.outer_product(weight_grad) * deformationGradientElastic;

    // // Check to make sure we should do these calculations?
    // if (del_elastic[0][0] < MATRIX_EPSILON &&
    //     del_elastic[0][1] < MATRIX_EPSILON &&
    //     del_elastic[1][0] < MATRIX_EPSILON &&
    //     del_elastic[1][1] < MATRIX_EPSILON)
    //     return Vector2f(0);

    // // Compute R^T*dF - dF^TR
    // // It is skew symmetric, so we only need to compute one value (three for
    // 3D) float y =
    //     (polar_r[0][0] * del_elastic[1][0] +
    //      polar_r[1][0] * del_elastic[1][1]) -
    //     (polar_r[0][1] * del_elastic[0][0] + polar_r[1][1] *
    //     del_elastic[0][1]);
    // // Next we need to compute MS + SM, where S is the hermitian matrix
    // // (symmetric for real valued matrices) of the polar decomposition and M
    // is
    // // (R^T*dR); This is equal to the matrix we just found (R^T*dF ...), so
    // we
    // // set them equal to eachother Since everything is so symmetric, we get a
    // // nice system of linear equations once we m->multiply everything out.
    // (see
    // // pdf for details) In the case of 2D, we only need to solve for one
    // // variable (three for 3D)
    // float x = y / (polar_s[0][0] + polar_s[1][1]);
    // // Final computation is deltaR = R*(R^T*dR)
    // Matrix2f del_rotate = Matrix2f(-polar_r[1][0] * x, polar_r[0][0] * x,
    //                                -polar_r[1][1] * x, polar_r[0][1] * x);

    // // We need the cofactor matrix of F, JF^-T
    // Matrix2f cofactor = deformationGradientElastic.cofactor();

    // // The last matrix we need is delta(JF^-T)
    // // Instead of doing the complicated matrix derivative described in the
    // paper
    // // we can just take the derivative of each individual entry in JF^-T;
    // JF^-T
    // // is the cofactor matrix of F, so we'll just hardcode the whole thing
    // For
    // // example, entry [0][0] for a 3x3 matrix is 	cofactor = e*i - f*h
    // // derivative
    // //= (e*Di + De*i) - (f*Dh + Df*h) 	where the derivatives (capital D)
    // come
    // // from our precomputed delta(F) In the case of 2D, this turns out to be
    // // just
    // // the cofactor of delta(F) For 3D, it will not be so simple
    // Matrix2f del_cofactor = del_elastic.cofactor();

    // // Calculate "A" as given by the paper
    // // Co-rotational term
    // Matrix2f Ap = del_elastic - del_rotate;
    // Ap *= 2 * m->mu;
    // // Primary contour term
    // cofactor *= cofactor.frobeniusInnerProduct(del_elastic);
    // del_cofactor *= (deformationGradientElastic.determinant() - 1);
    // cofactor += del_cofactor;
    // cofactor *= m->lambda;
    // Ap += cofactor;

    // // Put it all together
    // // Parentheses are important; M*M*V is slower than M*(M*V)
    // return volume *
    //        (Ap * (deformationGradientElastic.transpose() * weight_grad));
}

// #endif

SnowParticleSet::SnowParticleSet() : particles()
{
}

SnowParticleSet::~SnowParticleSet()
{
    // TODO check if this destroy way is okay
    for (auto& oneParticle : particles)
    {
        delete (oneParticle);
    }
}

void SnowParticleSet::addParticle(SnowParticle* sp)
{
    particles.push_back(sp);
}

void SnowParticleSet::addParticle(const Vector3f& pos, SnowParticleMaterial* m)
{
    SnowParticle* sp = new SnowParticle(pos, m);
    particles.push_back(sp);
}

void SnowParticleSet::addParticlesInAShape(Shape* s, SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    // std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, m);
        }
    }
}

void SnowParticleSet::update()
{
    maxVelocity = 0;
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i]->updatePos();
        particles[i]->updatePureElasticGradient();
        particles[i]->updateCombinedPElasticGradient();
        // Update max velocity, if needed
        float vel = particles[i]->velocity.norm();
        if (vel > maxVelocity) maxVelocity = vel;
    }
}