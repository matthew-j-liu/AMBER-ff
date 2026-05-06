#include <gtest/gtest.h>

#include "read_xyz.hpp"
#include <stdexcept>
#include <cmath>


// ASSERTIONS (weird input files)

// Missing Input Coord File
TEST(read_xyz, ThrowsMissingFile) {
    EXPECT_THROW(read_xyz("input_molecules/missing_file.xyz"), std::runtime_error);
}

// Atom count does not match num of atom rows
TEST(read_xyz, ThrowsBadAtomCount) {
    EXPECT_THROW(read_xyz(), std::runtime_error);
}

// non numerical coordinates
TEST(read_xyz, ThrowsBadAtomCount) {
    EXPECT_THROW(read_xyz(), std::runtime_error);
}


// Atom not an element in list of acceptable elements 
TEST(read_xyz, ThrowsBadAtomCount) {
    EXPECT_THROW(read_xyz(), std::runtime_error);
}


// ZERO ATOMS
// create a new ip file for this test case 
TEST(read_xyz, HasOneAtomOnly) {
    auto mol = read_xyz("h_atom.xyz");  // add single h atom file to input molecules 
    EXPECT_EQ(mol.num_atoms(), 0);
    EXPECT_EQ(mol.num_b 

// create a new ip file for this test case
TEST(read_xyz, HasOneAtomOnly) {
    auto mol = read_xyz("h_atom.xyz");  // add single h atom file to input molecules 
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.num_bonds(), 0);
    EXPECT_EQ(mol.atoms[0].element, "C");
    EXPECT_EQ(mol.atoms[0].position[0], 0.0);
    EXPECT_EQ(mol.atoms[0].position[1], 0.0);
    EXPECT_EQ(mol.atoms[0].position[2], 0.0);
}

TEST(read_xyz, HasAtleastOneAtom) {
    auto mol = read_xyz(MOL + "ethane.xyz");
}

/* Ethane and Beyond
1. number of atoms
2. number of bonds 
3. Elements ID
4. Positions 
5. Amber Type
6. Neighbors in graph (or should this be in MoleculeGraph struct test cases??)
*/ 

TEST(read_xyz, IsEthaneNumAtoms) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    EXPECT_EQ(mol.num_bonds(), 7);
}

TEST(read_xyz, IsEthaneNumBonds) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    EXPECT_EQ(mol.num_bonds(), 7);
}

TEST(read_xyz, IsEthaneElementID) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    EXPECT_EQ(mol.atoms[0].element, "C");
    EXPECT_EQ(mol.atoms[1].element, "C");
    for (int i = 2; i < 8; i++)
        EXPECT_EQ(mol.atoms[i].element, "H");
}

TEST(read_xyz, EthaneCarbonPosition) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    EXPECT_EQ(mol.atoms[0].position[0],  0.0);
    EXPECT_EQ(mol.atoms[0].position[1],  0.761550);
    EXPECT_EQ(mol.atoms[0].position[2],  0.0);
}

TEST(read_xyz, IsEthaneAmberTypeAssignment) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    EXPECT_EQ(mol.atoms[0].amber_type, "c3");
    EXPECT_EQ(mol.atoms[1].amber_type, "c3");
    for (int i = 2; i < 8; i++) 
        EXPECT_EQ(mol.atoms[i].amber_type, "hc");
}

TEST(read_xyz, EthaneCheckNeighbors) {
    auto mol = read_xyz(MOL + "ethane.xyz");
    //// 
}

