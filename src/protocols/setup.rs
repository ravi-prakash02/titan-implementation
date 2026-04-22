use rand::prelude::StdRng;
use rand::SeedableRng;
use crate::group_whir_committer::{GroupWhirProverMerkleState, GroupWhirVerifierMerkleState, MerkleProofStrategy};
use crate::merkle_tree::{Sha256Compress, Sha256LeafHash, Sha256MerkleTreeParams};
use crate::protocols::groupbulletproof::BulletProofParams;
use crate::protocols::titanpcs::{TitanPolyCommitment, TitanProverParams, TitanSetup, TitanSetupConfig, TitanVerifierParams};
use crate::utils::create_smooth_domain;

type PallasG = crate::pastatypes::Point;
type PallasF = crate::pastatypes::Scalar;

type ArkPallasG = crate::arkpallastypes::ArkPallasPoint;


impl TitanPolyCommitment<PallasG, PallasF> {
    pub fn setup(config: &TitanSetupConfig) -> (TitanProverParams<PallasG, PallasF>, TitanVerifierParams<PallasG, PallasF>) {
        let mut rng = StdRng::from_entropy();
        let domain_g1 = create_smooth_domain(config.domain_g1_size);
        let domain_g2 = create_smooth_domain(config.domain_g2_size);
        let bp_setup = BulletProofParams::<PallasG, PallasF>::new(config.m - config.m1, config.l2);
        let titan_setup = TitanSetup {
            m: config.m,
            m1: config.m1,
            l1: config.l1,
            num_queries: config.num_queries,
            domain_g1: domain_g1.1.clone(),
            bp_setup: bp_setup,
        };



        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<PallasG, Sha256LeafHash<PallasG>, Sha256Compress>(&mut rng);
        //let (leaf_hash_params_prover, two_to_one_params_prover) =
        //    crate::merkle_tree::default_config::<crate::pastatypes::Point, Sha256LeafHash<crate::pastatypes::Point>, Sha256Compress>(&mut rng);


        let verifier_state = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<PallasG>, PallasG>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let verifier_state_prover = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<PallasG>, PallasG>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let verifier_params = TitanVerifierParams {
            pp: titan_setup.clone(),
            state: verifier_state
        };

        let prover_params = TitanProverParams {
            pp: titan_setup.clone(),
            state: GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed),
            verifier_state: verifier_state_prover,
        };

        (prover_params, verifier_params)
    }
}

impl TitanPolyCommitment<ArkPallasG, PallasF> {
    pub fn setup(config: &TitanSetupConfig) -> (TitanProverParams<ArkPallasG, PallasF>, TitanVerifierParams<ArkPallasG, PallasF>) {
        let mut rng = StdRng::from_entropy();
        let domain_g1 = create_smooth_domain(config.domain_g1_size);
        let domain_g2 = create_smooth_domain(config.domain_g2_size);
        let bp_setup = BulletProofParams::<ArkPallasG, PallasF>::new(config.m - config.m1, config.l2);
        let titan_setup = TitanSetup {
            m: config.m,
            m1: config.m1,
            l1: config.l1,
            num_queries: config.num_queries,
            domain_g1: domain_g1.1.clone(),
            bp_setup: bp_setup,
        };



        let (leaf_hash_params, two_to_one_params) =
            crate::merkle_tree::default_config::<crate::arkpallastypes::ArkPallasPoint, Sha256LeafHash<crate::arkpallastypes::ArkPallasPoint>, Sha256Compress>(&mut rng);
        //let (leaf_hash_params_prover, two_to_one_params_prover) =
        //    crate::merkle_tree::default_config::<crate::pastatypes::Point, Sha256LeafHash<crate::pastatypes::Point>, Sha256Compress>(&mut rng);


        let verifier_state = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<crate::arkpallastypes::ArkPallasPoint>, crate::arkpallastypes::ArkPallasPoint>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let verifier_state_prover = GroupWhirVerifierMerkleState::<Sha256MerkleTreeParams<crate::arkpallastypes::ArkPallasPoint>, crate::arkpallastypes::ArkPallasPoint>::new(
            MerkleProofStrategy::Compressed,
            leaf_hash_params,
            two_to_one_params,
        );

        let verifier_params = TitanVerifierParams {
            pp: titan_setup.clone(),
            state: verifier_state
        };

        let prover_params = TitanProverParams {
            pp: titan_setup.clone(),
            state: GroupWhirProverMerkleState::new(MerkleProofStrategy::Compressed),
            verifier_state: verifier_state_prover,
        };

        (prover_params, verifier_params)
    }
}
